classdef ShapePredictor
    %SHAPEPREDICTOR Predict shape from measurement
    
    properties (Constant)
        EPSILON = 10^(-12);
        DEFAULT_PVEC0 = [0;0;0]
        DEFAULT_tvec0 = [1;0;0]
        DEFAULT_bvec0 = {[0;1;0],[0;0;1]}
    end
    
    methods (Static)
        function lp_data = lowpassFilter(data, beta, fc)
            %LOWPASSFILTER Apply kaiser window function and lowpass filter
            %{
            Args:
                data (Numeric): Data.
                beta (Numeric): Parameter of Kaiser window function.
                fc (Numeric): Normalized passband frequency.
            Returns:
                lp_data (Numeric): Data with lowpass fileter applied.
            %}
            
            L = length(data);
            w = kaiser(L, beta);
            offset = mean(data,2);
            x = transpose(data - offset).*w;
            [y, ~] = lowpass(x, fc,'Steepness',0.95,'StopbandAttenuation',30);
            lp_data = transpose(y./w) + offset;
        end
        
        function meanAngle = meanAngle(angles,dim)
            %MEANANGLE Get mean angle
            %{
            Args:
                angles (Numeric): Angles.
                dim (Numeric): Dimensiion to average.
            Returns:
                meanAngle (Numeric): Mean angle.
            %}
            
            meanAngle = atan2(sum(sin(angles),dim),sum(cos(angles),dim));
        end
        
        function maAngle = movmeanAngle(angles,dim,windowLen)
            %MOVMEANANGLE Apply Moving Average to angles
            %{
            Args:
                angles (Numeric): Angles.
                dim (Numeric): Dimensiion to average.
                windowLen (Numeric): Sample num of window for Moving Average.
            Output:
                maAngle (Numeric): Angle with Moving Average applied.
            %}
            
            if dim == 2
                angles = transpose(angles);
            end
            maAngle = zeros(size(angles));
            nu = ceil((windowLen - 1)/2);
            nl = floor((windowLen - 1)/2);
            [s1,s2] = size(angles);
            for i = 1:s1
                maAngle(i,:) = mcfShapePredictor.ShapePredictor.meanAngle(angles(max([1 i-nu]):min([s2 i+nl]),:),1);
            end
            if dim == 2
                maAngle = transpose(maAngle);
            end
        end
        
        function [s,pvec,prediction] = predict(result,varargin)
            % PREDICT
            %{
            Args:
                mcfParameters (Struct): Meta information about MCF.
                    PoissonRatio (Numeric): Poisson Ratio of MCF.
                    distance (Numeric): Distances between center and each core. Units: [micro m]
                    initialAngle (Numeric): Initial angls of each core. Units: [rad]
                    pitch (Function hundle): Pitches of each core. Units: [1/m]
                measurement (Struct): Measurement data.
                    s (Numeric): Arc length. Units: [m]
                    strain (Numeric): Measurement strain of deformed shape. Units: [strain]
                    initial_pvec (Numeric): Intial pvec
                    initial_tvec (Numeric): Intial tvec
                residualSrain (Numeric, Optional): Residual strain. Units: [strain]
                    This equals to measurement strain of straight shape without tension and torsion
                    , or difference true strain and measurement strain of any shape known.
                lpf (Empty or Struct, Parameter): If not empty, lowpass filter is apllied. This has follwing fields.
                    beta (Numeric): Parameter of Kaiser window function.
                    fc (Numeric): Normalized passband frequency.
                noTorsion (Logical, Parameter): If true, skip torsion calculation and set torsion strain and specificangle of twist zero.
                usedCores (Numeric, Parameter): The number of cores used to predict shape.
            Returns:
                s (Numeric): Arc length on MCF, or arc length. Units: [m]
                pvec (Numeric): Position vector.
                prediction (Struct): Prediction of prediction process.
            %}
            
            % Options
            default_residualStrain = [];
            default_lpf = [];
            default_noTorsion = false;
            default_usedCores = [];
            
            p = inputParser;
            addOptional(p, 'residualStrain', default_residualStrain);
            addParameter(p, 'lpf', default_lpf);
            addParameter(p, 'noTorsion', default_noTorsion);
            addParameter(p, 'usedCores', default_usedCores);
            parse(p, varargin{:});
            
            prediction = struct();
            s = result.measurement.s;
            
            % Set strain
            assert(all(size(result.measurement.strain)==size(p.Results.residualStrain)),...
                'Size of residual strain does not match that of measurement strain.')
            [measurementStrain,residualStrain] = mcfShapePredictor.ShapePredictor.setStrain(...
                result.measurement.strain,p.Results.residualStrain,p.Results.lpf);
            measurementStrain = vertcat(measurementStrain,mean(measurementStrain));  %追加
            residualStrain = vertcat(residualStrain,mean(residualStrain)); %追加
            % Set properties
            
            %訂正あり
            [coreNum,sampleNum] = size(measurementStrain);
            coreNum =coreNum + 1;      %追加
            sampleNum = sampleNum+1;   %追加
            PoissonRatio = result.mcfParameters.PoissonRatio;
            distance = result.mcfParameters.distance * 10^(-6);
            distance = horzcat(distance,0);   %追加
            initialAngle = result.mcfParameters.initialAngle;
            initialAngle = horzcat(initialAngle,0); %追加
            
            if iscell(result.mcfParameters.pitch)
                pitch = zeros(coreNum-1,length(s));
%                 for i = 1:coreNum
%                     pitch_fun = result.mcfParameters.pitch{i};
%                     pitch(i,:) = pitch_fun(s);
%                 end
             
            else
                pitch = repmat(transpose(result.mcfParameters.pitch),1,length(s));
            end
            
            initialPvec = result.measurement.initialPvec;
            initialTvec = result.measurement.initialTvec;
            initialBeta = result.measurement.initialBeta;
            
            % usedCores option
%             if ~isempty(p.Results.usedCores)
%                 measurementStrain = measurementStrain(p.Results.usedCores,:);
%                 residualStrain = residualStrain(p.Results.usedCores,:);
%                 
%                 coreNum = length(p.Results.usedCores);
%                 distance = distance(p.Results.usedCores);
%                 initialAngle = initialAngle(p.Results.usedCores);
%                 pitch = pitch(p.Results.usedCores,:);
%             end

            %卒論用訂正版
            if ~isempty(p.Results.usedCores)
                measurementStrain = measurementStrain(p.Results.usedCores,:);
                measurementStrain = vertcat(measurementStrain,mean(measurementStrain));
                residualStrain = residualStrain(p.Results.usedCores,:);
                residualStrain = vertcat(residualStrain,mean(residualStrain));
                
                coreNum = length(p.Results.usedCores)+1;
                distance = distance(p.Results.usedCores);
                disp(distance);
                initialAngle = initialAngle(p.Results.usedCores);
                initialAngle = horzcat(initialAngle,0);
                pitch = pitch(p.Results.usedCores,:);
                
            end
            %訂正終了
            
            centerCore = distance == 0;   %コアの中心からの距離が0のものがあればtrue
            
            
            assert(any(centerCore),'No center core.')    %assertはfalseならエラーをスロー
            % assert(sum(centerCore)==1,'Too many center cores.')
            
            % Sprit strain
            %   Residual Strain
            strain_tmp = measurementStrain - residualStrain;  %残留歪み取り除き
            %   Axial Strain
            axialStrain = mcfShapePredictor.ShapePredictor.calAxialStrain(...
                strain_tmp,PoissonRatio,distance,pitch,centerCore);
            strain_tmp = strain_tmp - axialStrain;   %軸歪み取り除き
            % Check condition for classification
            condition = mcfShapePredictor.ShapePredictor.checkCondition(...
                p.Results.noTorsion,distance,pitch,initialAngle,centerCore);
            switch condition
                case {1,2}  % Torsion strain is splited easily.
                    %   Torsion Strain
                    if condition == 1  % No torsion or cannot mesure torsion
                        torsionStrain = zeros(coreNum-1,sampleNum-1);
                        saot = zeros(1,sampleNum);
                    else  % Bending strains cancel each other
                        [saot,torsionStrain] = mcfShapePredictor.ShapePredictor.calSaot(strain_tmp,distance,pitch,centerCore);
                    end
                    %   Bending Strain
                    
                    bendingStrain = strain_tmp - torsionStrain;
                    
                    % Calculate Crv and Beta
                    [crv,beta] = mcfShapePredictor.ShapePredictor.calCrvAndBeta(...
                        s,bendingStrain,PoissonRatio,initialBeta,distance,pitch,initialAngle,centerCore,saot);
                    
                case 3  % Use solver function.
                    % Cal Crv, Beta ans Saot with solver.
                    [crv,beta,saot,fvals] = mcfShapePredictor.ShapePredictor.calCrvBetaAndSaot(...
                        s,strain_tmp,PoissonRatio,distance,pitch,initialAngle,centerCore);
                    prediction.fvals = fvals;
                    %   Torsion and Bending Strain
                    [torsionStrain,bendingStrain] = mcfShapePredictor.ShapePredictor.calTorsionAndBendingStrain(...
                        s,crv,beta,saot,PoissonRatio,distance,pitch,initialAngle);
                    
                otherwise
                    error('Undifined condition.')
            end
            [ds,trs] = mcfShapePredictor.ShapePredictor.calDsAndTrs(s,beta);
            
            prediction.s = s;
            prediction.strain = measurementStrain;
            prediction.residualStrain = residualStrain;
            prediction.axialStrain = axialStrain;
            prediction.torsionStrain = torsionStrain;
            prediction.bendingStrain = bendingStrain;
            prediction.saot = saot;
            prediction.crv = crv;
            prediction.beta = beta;
            prediction.ds = ds;
            prediction.trs = trs;
            
            [pvec,tvec,bvec,nvec] = mcfShapePredictor.ShapePredictor.coords(ds,crv,beta,trs,initialPvec,initialTvec);
            
            prediction.pvec = pvec;
            prediction.tvec = tvec;
            prediction.bvec = bvec;
            prediction.nvec = nvec;
        end
        
        function [measurementStrain,residualStrain] = setStrain(measurementStrain,residualStrain,lpf)
            %SETSTRAIN Set strain and apply lowpass filter.
            %{
            Args:
                measurementStrain (Numeric): Measurement strain of deformed shape. Units: [strain]
                residualSrain (Numeric): Residual strain. Units: [strain]
                    This equals to measurement strain of straight shape without tension and torsion
                    , or difference true strain and measurement strain of any shape known.
                lpf (Empty or Struct): If empty, lowpass filter is not apllied. This has follwing fields.
                    beta (Numeric): Parameter of Kaiser window function.
                    fc (Numeric): Normalized passband frequency.
            Returns:
                measurementStrain (Numeric): Measurement strain of deformed shape with lowpass filter allpied. Units: [strain]
                residualSrain (Numeric): Residual strain with lowpass filter allpied. Units: [strain]
            %}
            
            if ~isempty(lpf)
                measurementStrain = mcfShapePredictor.ShapePredictor.lowpassFilter(measurementStrain,lpf.beta,lpf.fc);
            end
            if isempty(residualStrain)
                residualStrain = zeros(size(measurementStrain));
            elseif ~isempty(lpf)
                residualStrain = mcfShapePredictor.ShapePredictor.lowpassFilter(residualStrain,lpf.beta,lpf.fc);
            end
        end
        
        function axialStrain = calAxialStrain(strain_atb,PoissonRatio,distance,pitch,centerCore)
            %CALAXIALSTRAIN Calculate axial strain
            %{
            Args:
                strain_atb (Numeric): Sum of axial strain, torsion strain, and bending strain. Units : [strain]
                PoissonRatio (Numeric): Poisson ratio of MCF.
                distance (Numeric): Distances between center and each core. Units: [micro m]
                pitch (Numeric): Pitches of each core. Units: [1/m]
                centerCore (Logical): Index of center core.
            Returns:
                axialStrain (Numeric): Axial strain.
            %}
            
            coreNum = size(strain_atb,1);
           
            k1 = (1 - PoissonRatio.*(2*pi*pitch.*transpose(distance)))./hypot(2*pi*pitch.*transpose(distance),1);
            axialStrain = k1.*repmat(strain_atb(centerCore,:), coreNum, 1);
        end
        
        function condition = checkCondition(noTorsion,distance,pitch,initialAngle,centerCore)
            %CHECKCONDITION Check condition for classification
            %{
            Args:
                noTorsion (Logical): If true, skip torsion calculation and set torsion strain and specificangle of twist zero.
                distance (Numeric): Distances between center and each core. Units: [micro m]
                pitch (Numeric): Pitches of each core. Units: [1/m]
                initialAngle (Numeric): Initial angls of each core. Units: [rad]
                centerCore (Logical): Index of center core.
            Returns:
                condition (Numeric): Condition of parameters.
                    Condition 1: No torsion or cannot mesure torsion.
                    Condition 2: Bending strains cancel each other.
                    Condition 3: Use solver function to split torsion and bending strain.
            %}
            
            if noTorsion || all(distance==0) || all(pitch==0,'all')
                condition = 1;
                fprintf('Condition 1: No torsion or cannot mesure torsion.\n')
            else
                const_distance = length(unique(distance(~centerCore))) == 1;
                const_pitch = length(unique(pitch(~centerCore,:))) == 1;
                
                initialAngle_nc = initialAngle(~centerCore);
                dia = initialAngle_nc - initialAngle_nc(1);
                appropriateInitialAngle = abs(sum(sin(dia))) < mcfShapePredictor.ShapePredictor.EPSILON...
                    && abs(sum(cos(dia))) < mcfShapePredictor.ShapePredictor.EPSILON;
                
                if (const_distance && const_pitch) && appropriateInitialAngle
                    condition = 2;
                    fprintf('Condition 2: Bending strains cancel each other.\n')
                else
                    condition = 3;
                    warning('Condition 3: Use solver function to split torsion and bending strain.')
                end
            end
        end
        
        function [saot,torsionStrain] = calSaot(strain_tb,distance,pitch,centerCore)
            %CALSAOT Calculate Specific angle of twist and torsion strain
            %{
            Args:
                strain_tb (Numeic): Sum of torsion strain and bending strain. Units: [strain]
                distance (Numeric): Distances between center and each core. Units: [micro m]
                pitch (Numeric): Pitches of each core. Units: [1/m]
                centerCore (Logical): Index of center core.
            Returns:
                saot (Numeric): Specific angle of twist. Units: [rad/m]
                torsionStrain (Numeric): Torsion strian. Units: [strain]
            %}
            
            [coreNum, sampleNum] = size(strain_tb);
            k2 = (2*pi.*pitch.*transpose(distance))./hypot(2*pi.*pitch.*transpose(distance),1);
            torsionStrain = zeros(coreNum,sampleNum);
            coef = k2.*transpose(distance);
            
            strain_tb_nc = strain_tb(~centerCore,:);
            torsionStrain(~centerCore,:) = repmat(mean(strain_tb_nc),sum(~centerCore),1);  % Cancel bending strains
            saot = mean(strain_tb_nc)./mean(coef(~centerCore),'all');
        end
        
        function domega = domega(s,initialAngle,pitch,core1,core2)
            %DOMEGA Get difference of omega between core1 and core2
            %{
            Args:
                s (Numeric): Arc length. Units: [m]
                initialAngle (Numeric): Initial angls of each core. Units: [rad]
                pitch (Numeric): Pitches of each core. Units: [1/m]
                core1 (Numeric): Index of core 1.
                core2 (Numeric): Index of core 2.
            Returns:
                domega (Numeric): Difference of omega between core1 and core2.
            %}
            
            domega = initialAngle(core1) - initialAngle(core2) + 2*pi*s.*(pitch(core1,:) - pitch(core2,:));
        end
        
        function omega = omega(s,initialAngle,pitch,saot)
            %OMEGA  Get omega, angle of core at the arc length.
            %{
            Args:
                s (Numeric): Arc length. Units: [m]
                initialAngle (Numeric): Initial angls of each core. Units: [rad]
                pitch (Numeric): Pitches of each core. Units: [1/m]
                saot (Numeric): Specific angle of twist. Units: [rad/m]
            Returns:
                omega (Numeric): Omega, angle of core at the sition.
            %}
            
            torsionAngle = [0 cumsum(saot(1,1:end-1).*diff(s))];
            disp(torsionAngle);
            omega = transpose(initialAngle) + 2*pi*s.*pitch + torsionAngle;
        end
        
        function [crv,beta] = calCrvAndBeta(s,bendingStrain,PoissonRatio,...
                initialBeta,distance,pitch,initialAngle,centerCore,saot)
            %CALCRVANDBETA Calculate Curvature and Beta
            %{
            Args:
                s (Numeric): Arc length. Units: [m]
                bendingStrain (Numeric): Bending strain. Units: [strain]
                PoissonRatio (Numeric): Poisson ratio of MCF.
                initialBeta (Numeric): Initial Beta.
                distance (Numeric): Distances between center and each core. Units: [micro m]
                pitch (Numeric): Pitches of each core. Units: [1/m]
                initialAngle (Numeric): Initial angls of each core. Units: [rad]
                centerCore (Logical): Index of center core.
                saot (Numeric): Specific angle of twist. Units: [rad/m]
            Returns:
                crv (Numeric): Curvature. Units: [/m]
                beta (Numeric): Beta, or benging angle. Units: [rad]
            %}
            
            % Remove center core
            bendingStrain = bendingStrain(~centerCore,:);
            distance = distance(~centerCore);
            initialAngle = initialAngle(~centerCore);
            pitch = pitch(~centerCore,:);
            [~,sampleNum] = size(bendingStrain);
            
            k1 = (1 - PoissonRatio.*(2*pi*pitch.*transpose(distance)))./hypot(2*pi*pitch.*transpose(distance),1);
            omega = mcfShapePredictor.ShapePredictor.omega(s,initialAngle,pitch,saot);
            A = cat(3,transpose(k1.*transpose(distance).*cos(omega)),transpose(k1.*transpose(distance).*sin(omega)));
            
            x = zeros(sampleNum,2);
            if size(bendingStrain,1)==2
                for i = 1:sampleNum
                    x(i,:) = squeeze(A(i,:,:))\squeeze(bendingStrain(:,i));
                end
            else
                for i = 1:sampleNum
                    A_ = squeeze(A(i,:,:));
                    e_ = squeeze(bendingStrain(:,i));
                    x(i,:) = (transpose(A_)*A_)\transpose(A_)*e_;
                end
            end
            crv = vecnorm(transpose(x));
            bending_idx = crv ~= 0;
            % beta in front of bending
            beta = ones(1,sampleNum)*initialBeta;
            % beta in bending
            beta(bending_idx) = atan2(x(bending_idx,2),x(bending_idx,1));
            % beta behind bending
            if any(~bending_idx)
                dist = s(~bending_idx) - transpose(s(bending_idx));
                dist(dist <= 0) = nan;
                [dist, idx_] = min(dist,[],1);
                bending_idx_ = find(bending_idx);
                not_bending_idx = find(~bending_idx);
                beta(not_bending_idx(~isnan(dist))) = beta(bending_idx_(idx_(~isnan(dist))));
            end            
        end
        
        function [crv,beta,saot,fvals] = calCrvBetaAndSaot(s,strain_tb,PoissonRatio,distance,pitch,initialAngle,centerCore)
            %CALCRVBEATAANDSAOT Calculate curvature, beta, and specific angle of twist with solver function
            %{
            Solve following equations about crv, beta, and saot.
                    e_i = k1_i * r_i * crv * cos(omega - beta) + k2_i * r * saot
                    omega = a_i + Integrate 2 * pi * pitch_i + saot(t) dt from 0 to s.
            Args:
                s (Numeric): Arc length. Units: [m]
                strain_tb (Numeric): Sum of torsion strain and bending strain. Units: [strain]
                PoissonRatio (Numeric): Poisson ratio of MCF.
                distance (Numeric): Distances between center and each core. Units: [micro m]
                pitch (Numeric): Pitches of each core. Units: [1/m]
                initialAngle (Numeric): Initial angls of each core. Units: [rad]
                centerCore (Logical): Index of center core.
            Returns:
                crv (Numeric): Curvature. Units: [/m]
                beta (Numeric): Beta, or benging angle. Units: [rad]
                saot (Numeric): Specific angle of twist. Units: [rad/m]
                fvals (Numeric): Residual erros of optimization. Units: [strain]
            %}
            
            sampleNum = size(strain_tb,2);
            crv = zeros(1,sampleNum);
            beta = zeros(1,sampleNum);
            saot = zeros(1,sampleNum);
            fvals = zeros(sum(~centerCore),sampleNum);
            
            % Remove center core
            strain_tb = strain_tb(~centerCore,:);
            distance = distance(~centerCore);
            initialAngle = initialAngle(~centerCore);
            pitch = pitch(~centerCore,:);
            
            k1 = (1 - PoissonRatio.*(2*pi*pitch.*transpose(distance)))./hypot(2*pi*pitch.*transpose(distance),1);
            k2 = (2*pi.*pitch.*transpose(distance))./hypot(2*pi*pitch.*transpose(distance),1);
            
            function res = objectiveFun(crv,beta,saot,saotHist,s,strain_tb,k1,k2,distance,initialAngle,pitch)
                %OBJECTIVEFUN Objective function
                
                omega = mcfShapePredictor.ShapePredictor.omega(s,initialAngle,pitch,[saotHist,saot]);
                bendingStrain = k1.*transpose(distance)*crv.*cos(omega(:,end) - beta);
                torsionStrain = k2.*transpose(distance)*saot;
                res = bendingStrain + torsionStrain - strain_tb;
            end
            
            x = zeros(1,3);
            for i = 1:sampleNum
                fun = @(x) objectiveFun(x(1),x(2),x(3),saot(1:i-1),s(1:i),strain_tb(:,i),k1(:,i),k2(:,i),distance,initialAngle,pitch(:,i));
                x0 = x;
                options = optimoptions('fsolve',...
                    'Display','off',...
                    'OptimalityTolerance',mcfShapePredictor.ShapePredictor.EPSILON,...
                    'FunctionTolerance',mcfShapePredictor.ShapePredictor.EPSILON);
                [x,fval,exitflag,output] = fsolve(fun,x0,options);
                crv(i) = x(1);
                beta(i) = mcfShapePredictor.Utils.unwrap(x(2));
                saot(i) = x(3);
                fvals(:,i) = fval;
            end
        end
        
        function [torsionStrain,bendingStrain] = calTorsionAndBendingStrain(...
                s,crv,beta,saot,PoissonRatio,distance,pitch,initialAngle)
            %CALTORSIONANDBENDINGSTRAIN Calculate Torsion and Bending strain
            %{
            Args:
                s (Numeric): Arc length. Units: [m]
                crv (Numeric): Curvature. Units: [/m]
                beta (Numeric): Beta, or benging angle. Units: [rad]
                saot (Numeric): Specific angle of twist. Units: [rad/m]
                PoissonRatio (Numeric): Poisson ratio of MCF.
                distance (Numeric): Distances between center and each core. Units: [micro m]
                pitch (Numeric): Pitches of each core. Units: [1/m]
                initialAngle (Numeric): Initial angls of each core. Units: [rad]
                centerCore (Logical): Index of center core.
            Returns:
                torsionStrain (Numeric): Torsion strian. Units: [strain]
                bendingStrain (Numeric): Bending strian. Units: [strain]
            %}
            
            k1 = (1 - PoissonRatio.*(2*pi*pitch.*transpose(distance)))./hypot(2*pi*pitch.*transpose(distance),1);
            k2 = (2*pi.*pitch.*transpose(distance))./hypot(2*pi*pitch.*transpose(distance),1);
            omega = mcfShapePredictor.ShapePredictor.omega(s,initialAngle,pitch,saot);
            
            torsionStrain = k2.*transpose(distance).*saot;
            bendingStrain = k1.*transpose(distance).*crv.*cos(omega - beta);
        end
        
        function [ds,trs] = calDsAndTrs(s,beta)
            %CALDSANDTRS Calculate ds and trs
            %{
            Args:
                s (Numeric): Arc Length on MCF. Units: [m]
                beta (Numeric): Beta, or benging angle. Units: [rad]
            Returns:
                ds (Numeric): Derivative of arc length s.
                trs (Numeric): Derivative of beta.
            %}
            
            ds = [diff(s) 0];
            trs = [diff(beta)./diff(s) 0];
        end
        
        function [pvec,tvec,nvec,bvec] = coords(ds,crv,beta,trs,pvec0,tvec0)
            %COORDS Calculate coords
            %{
            Args:
                ds (Numeric): Derivative of arc length s.
                crv (Numeric): Curvature. Units: [/m]
                beta (Numetic): Beta, or benging angle. Units: [rad]
                trs (Numeric): Derivative of beta.
                pvec0 (Numeric): Initial pvec.
                tvec0 (Numeric): Initial tvec.
            Returns:
                pvec (Numeric): Position vec.
                tvec (Numeric): Tangent vec.
                nvec (Numeric): Normal vec.
                bvec (Numeric): Binormal vec.
            %}
            
            M = mcfShapePredictor.ShapePredictor.frenet_serret_matrix(crv,trs,ds);
            
            [nvec0,bvec0] = mcfShapePredictor.Utils.calNvec0AndBvec0(tvec0,beta(1));
            [pvec,tvec,nvec,bvec] = mcfShapePredictor.ShapePredictor.recur_frenet_serret(...
                M,pvec0,tvec0,nvec0,bvec0);
        end
        
        function M = frenet_serret_matrix(crv,trs,ds)
            %FRENET_SERRET_MATRIX frenet_serret matrix formulation
            %{
            Args:
                crv (Numeric): Curvature. Units: [/m]
                trs (Numeric): Curvature. Units: [rad/m]
                ds (Numeri): Difference of arg length s.
            Returns:
                M (Numeric): Frenet_serret matrix.
            %}
            
            theta = hypot(crv,trs).*ds;
            idx = crv ~= 0 | trs ~= 0;
            
            M = zeros(length(theta),4,4);
            M(:,1,1) = 1;
            M(idx,1,2) = crv(idx).^2.*sin(theta(idx))./hypot(crv(idx),trs(idx)) + trs(idx).^2.*ds(idx);
            M(idx,1,3) = crv(idx).*(1 - cos(theta(idx)));
            M(idx,1,4) = crv(idx).*trs(idx).*(ds(idx) - sin(theta(idx))./hypot(crv(idx),trs(idx)));
            M(idx,1,2:4) = M(idx,1,2:4)./reshape(hypot(crv(idx),trs(idx)).^2,[],1,1);
            M(~idx,1,2) = ds(~idx);
            M(idx,2,2) = crv(idx).^2.*cos(theta(idx)) + trs(idx).^2;
            M(idx,2,3) =  crv(idx).*hypot(crv(idx),trs(idx)).*sin(theta(idx));
            M(idx,2,4) = crv(idx).*trs(idx).*(1 - cos(theta(idx)));
            M(idx,3,2) = -crv(idx).*hypot(crv(idx),trs(idx)).*sin(theta(idx));
            M(idx,3,3) = (crv(idx).^2 + trs(idx).^2).*cos(theta(idx));
            M(idx,3,4) = trs(idx).*hypot(crv(idx),trs(idx)).*sin(theta(idx));
            M(idx,4,2) = crv(idx).*trs(idx).*(1 - cos(theta(idx)));
            M(idx,4,3) = -trs(idx).*hypot(crv(idx),trs(idx)).*sin(theta(idx));
            M(idx,4,4) = crv(idx).^2 + trs(idx).^2.*cos(theta(idx));
            M(idx,2:4,2:4) = M(idx,2:4,2:4)./reshape(hypot(crv(idx),trs(idx)).^2,[],1,1);
            M(~idx,2:4,2:4) = repmat(reshape(diag(ones(1,3)),1,3,3),sum(~idx),1,1);
        end
        
        function v = frenet_serret(M,v0)
            %FRENET_SERRET Frenet serret
            %{
            Args:
                M (Numeric): Frenet_serret matrix.
                v0 (Numeric): Initial vec.
            Returns
                v (Numeric): next vec.
            %}
            
            v(1:3:12) = M*v0(1:3:12);
            v(2:3:12) = M*v0(2:3:12);
            v(3:3:12) = M*v0(3:3:12);
            v(4:6) = v(4:6)/vecnorm(v(4:6));
            v(7:9) = v(7:9)/vecnorm(v(7:9));
            v(10:12) = v(10:12)/vecnorm(v(10:12));
        end
        
        function [pvec,tvec,nvec,bvec] = recur_frenet_serret(M,pvec0,tvec0,nvec0,bvec0)
            %RECUR_FRENET_SERRET Recuir frenet serret
            %{
            Args:
                M (Numeric): Frenet_serret matrix.
                pvec0 (Numeric): Initial pvec.
                tvec0 (Numeric): Initial tvec.
                nvec0 (Numeric): Initial nvec.
                bvec0 (Numeric): Initial bvec.
            Returns:
                pvec (Numeric): Position vec.
                tvec (Numeric): Tangent vec.
                nvec (Numeric): Normal vec.
                bvec (Numeric): Binormal vec.
            %}
            
            v = zeros(4*3,size(M,1));
            v(:,1) = [pvec0; tvec0; nvec0; bvec0];
            for i = 1:(size(M,1)-1)
                v(:,i+1) = mcfShapePredictor.ShapePredictor.frenet_serret(squeeze(M(i,:,:)),v(:,i));
            end
            pvec = v(1:3,:);
            tvec = v(4:6,:);
            nvec = v(7:9,:);
            bvec = v(10:12,:);
        end
        
        function save(filename,prediction)
            %SAVE Save prediction
            %{
            Args:
                filename (Char or String): File name.
                prediction (Struct): Prediction result.
            %}
            
            fields = fieldnames(prediction);
            for i = 1:length(fields)
                eval(strcat(fields{i}, ' = prediction.', fields{i}, ';'))
            end
            save(filename,fields{:})
        end
    end
end
