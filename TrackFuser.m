classdef TrackFuser < handle
    % TrackFuser 红外-雷达融合跟踪器类
    
    properties
        % 跟踪状态
        radarstate
        radarcov
        irstate
        ircov
        fusiontrace
        fusioncov
        is_initialized
    end
    
    methods
        function obj = TrackFuser()
            % TrackFuser 构造函数
            obj.reset();
        end
        
        function reset(obj)
            % reset 重置跟踪器状态
            obj.radarstate = zeros(6,1);
            obj.radarcov = zeros(6,6);
            obj.irstate = zeros(6,1);
            obj.ircov = zeros(6,6);
            obj.fusiontrace = zeros(6,1);
            obj.fusioncov = zeros(6,6);
            obj.is_initialized = true;
        end
        
        function track = update(obj, data)
            % update 更新跟踪状态
            % 输入: data - 1×5向量 [ir_az, ir_el, radar_az, radar_el, radar_r]
            % 输出: track - 结构体包含当前跟踪状态
            
            if ~obj.is_initialized
                obj.reset();
            end
            
            % 检查是否第一次调用（所有状态都为0）
            is_first_call = all(obj.radarstate(:,end) == 0) && ...
                           all(obj.irstate(:,end) == 0) && ...
                           all(obj.fusiontrace(:,end) == 0);
            
            % 提取当前状态
            state_r = obj.radarstate(:,end);
            if ndims(obj.radarcov) == 3
                cov_r = obj.radarcov(:,:,end);
            else
                cov_r = obj.radarcov;
            end
            
            state_ir = obj.irstate(:,end);
            if ndims(obj.ircov) == 3
                cov_ir = obj.ircov(:,:,end);
            else
                cov_ir = obj.ircov;
            end
            
            % 执行融合跟踪
            [state_r_update, cov_r_update, state_ir_update, cov_ir_update, fused_state, fused_cov] = ...
                obj.track_fusion(data, state_r, cov_r, state_ir, cov_ir, is_first_call);

            % 更新状态
            obj.radarstate = [obj.radarstate, state_r_update];
            obj.radarcov = cat(3, obj.radarcov, cov_r_update);
            obj.irstate = [obj.irstate, state_ir_update];
            obj.ircov = cat(3, obj.ircov, cov_ir_update);
            obj.fusiontrace = [obj.fusiontrace, fused_state];
            obj.fusioncov = cat(3, obj.fusioncov, fused_cov);

            % 返回当前跟踪状态
            track = struct(...
                'radarstate', obj.radarstate, ...
                'radarcov', obj.radarcov, ...
                'irstate', obj.irstate, ...
                'ircov', obj.ircov, ...
                'fusiontrace', obj.fusiontrace, ...
                'fusioncov', obj.fusioncov);
        end
    end
    
    methods (Static)
        function [state_r_update, cov_r_update, state_ir_update, cov_ir_update, fused_state, fused_cov] = ...
                track_fusion(data, state_r, cov_r, state_ir, cov_ir, is_first_call)
            
            %% ========== 统一门限 ==========
            reliabilityThreshold = chi2inv(0.9999, 9);
            epsilon = 1e-6;
            
            %% ========== 读取当前传感器数据 ==========
            ir_az = deg2rad(data(1));
            ir_el = deg2rad(data(2));
            radar_r = data(5);
            radar_az = deg2rad(data(3));
            radar_el = deg2rad(data(4));
            
            %% ========== 参数设置 ==========
            dt = 1;
            F = [1 0 0 dt 0  0;
                 0 1 0 0  dt 0;
                 0 0 1 0  0  dt;
                 0 0 0 1  0  0;
                 0 0 0 0  1  0;
                 0 0 0 0  0  1];
            Q = diag([1000, 50, 1000, 50, 1000, 50]);
            
            % 测量噪声
            sigma_r_radar = 1;
            sigma_theta_radar = deg2rad(1);
            sigma_theta_ir = deg2rad(0.2);
            
            R_r = diag([sigma_r_radar^2, sigma_theta_radar^2, sigma_theta_radar^2]);
            R_ir = diag([sigma_theta_ir^2, sigma_theta_ir^2]);
            
            %% ========== 初始化 ==========
            if is_first_call || all(state_r == 0)
                % 雷达初始化
                x_r = radar_r * cos(radar_el) * cos(radar_az);
                y_r = radar_r * cos(radar_el) * sin(radar_az);
                z_r = radar_r * sin(radar_el);
                state_r_update = [x_r; y_r; z_r; 0; 0; 0];
                cov_r_update = diag([100, 50, 100, 50, 100, 50]);
                
                % 红外初始化
                x_ir = radar_r * cos(ir_el) * cos(ir_az);
                y_ir = radar_r * cos(ir_el) * sin(ir_az);
                z_ir = radar_r * sin(ir_el);
                state_ir_update = [x_ir; y_ir; z_ir; 0; 0; 0];
                cov_ir_update = diag([100, 50, 100, 50, 100, 50]);
                
                % 初次融合
                fused_state = (state_r_update + state_ir_update) / 2;
                fused_cov = (cov_r_update + cov_ir_update) / 4;
                return;
            end
            
            %% ========== EKF预测 ==========
            [state_r_pred, cov_r_pred] = TrackFuser.ekf_predict(state_r, cov_r, F, Q);
            [state_ir_pred, cov_ir_pred] = TrackFuser.ekf_predict(state_ir, cov_ir, F, Q);
            
            %% ========== 干扰判定 ==========
            [zpred_r, H_r] = TrackFuser.predict_measurement_and_H(state_r_pred, 3);
            meas_r = [radar_r; radar_az; radar_el];
            
            y_r = meas_r - zpred_r;
            y_r(2) = mod(y_r(2) + pi, 2*pi) - pi;
            y_r(3) = mod(y_r(3) + pi, 2*pi) - pi;
            S_r = H_r * cov_r_pred * H_r' + R_r + epsilon * eye(size(R_r));
            mahal_r = y_r' / S_r * y_r;
            isRadarReliable = (mahal_r <= reliabilityThreshold);
            
            [zpred_ir, H_ir] = TrackFuser.predict_measurement_and_H(state_ir_pred, 2);
            meas_ir = [ir_az; ir_el];
            y_ir = meas_ir - zpred_ir;
            y_ir(1) = mod(y_ir(1) + pi, 2*pi) - pi;
            y_ir(2) = mod(y_ir(2) + pi, 2*pi) - pi;
            S_ir = H_ir * cov_ir_pred * H_ir' + R_ir + epsilon * eye(size(R_ir));
            mahal_ir = y_ir' / S_ir * y_ir;
            isIRReliable = (mahal_ir <= reliabilityThreshold);
            
            %% ========== EKF更新 ==========
            if isRadarReliable
                [state_r_update, cov_r_update] = TrackFuser.ekf_update(state_r_pred, cov_r_pred, meas_r, R_r);
            else
                state_r_update = state_r_pred;
                cov_r_update = cov_r_pred;
            end
            
            if isIRReliable
                [state_ir_update, cov_ir_update] = TrackFuser.ekf_update(state_ir_pred, cov_ir_pred, meas_ir, R_ir);
            else
                state_ir_update = state_ir_pred;
                cov_ir_update = cov_ir_pred;
            end
            
            %% ========== 状态融合 ==========
            [fused_state, fused_cov] = TrackFuser.fuse_tracks(...
                state_r_update, cov_r_update, state_ir_update, cov_ir_update, reliabilityThreshold);
        end
        
        function [z_pred, H] = predict_measurement_and_H(state, nm)
            x = state(1); y = state(2); z = state(3);
            r = sqrt(x^2 + y^2 + z^2);
            r_xy = sqrt(x^2 + y^2);
            
            if r == 0, r = eps; end
            if r_xy == 0, r_xy = eps; end
            
            if nm == 3
                z_pred = [r; atan2(y, x); atan2(z, r_xy)];
                H = zeros(3, 6);
                H(1,1) = x/r; H(1,2) = y/r; H(1,3) = z/r;
                H(2,1) = -y/(x^2+y^2); H(2,2) = x/(x^2+y^2);
                H(3,1) = -x*z/(r_xy*r^2); H(3,2) = -y*z/(r_xy*r^2); H(3,3) = r_xy/r^2;
                
            elseif nm == 2
                z_pred = [atan2(y, x); atan2(z, r_xy)];
                H = zeros(2, 6);
                H(1,1) = -y/(x^2+y^2); H(1,2) = x/(x^2+y^2);
                H(2,1) = -x*z/(r_xy*r^2); H(2,2) = -y*z/(r_xy*r^2); H(2,3) = r_xy/r^2;
            else
                error('nm必须为2或3');
            end
        end
        
        function [state_pred, cov_pred] = ekf_predict(state_prev, cov_prev, F, Q)
            state_pred = F * state_prev;
            cov_pred = F * cov_prev * F' + Q;
        end
        
        function [state_upd, cov_upd] = ekf_update(state_pred, cov_pred, measurement, R)
            nm = size(R,1);
            [z_pred, H] = TrackFuser.predict_measurement_and_H(state_pred, nm);
            
            y_res = measurement(:) - z_pred(:);
            
            if nm == 2
                y_res(1) = mod(y_res(1) + pi, 2*pi) - pi;
                y_res(2) = mod(y_res(2) + pi, 2*pi) - pi;
            else
                y_res(2) = mod(y_res(2) + pi, 2*pi) - pi;
                y_res(3) = mod(y_res(3) + pi, 2*pi) - pi;
            end
            
            S = H * cov_pred * H' + R;
            K = cov_pred * H' / S;
            state_upd = state_pred + K * y_res;
            I = eye(6);
            cov_upd = (I - K * H) * cov_pred;
        end
        
        function [fusedState, fusedCov] = fuse_tracks(stateRadar, covRadar, stateIR, covIR, reliabilityThreshold)
            posRadar = stateRadar(1:3);
            posIR = stateIR(1:3);
            covRadarPos = covRadar(1:3,1:3);
            covIRPos = covIR(1:3,1:3);
            
            diff = posRadar - posIR;
            maha = diff' / (covRadarPos + covIRPos) * diff;
            
            if maha < reliabilityThreshold
                W_r = inv(covRadar);
                W_i = inv(covIR);
                totalW = W_r + W_i;
                fusedState = totalW \ (W_r * stateRadar + W_i * stateIR);
                fusedCov = inv(totalW);
            else
                if trace(covRadar) < trace(covIR)
                    fusedState = stateRadar;
                    fusedCov = covRadar;
                else
                    fusedState = stateIR;
                    fusedCov = covIR;
                end
            end
        end
    end
end