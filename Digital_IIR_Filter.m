%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%                 Digital IIR Filter Creator                 %%%%
%%%%           Just run the program...Yes, it's that easy       %%%%
%%%%                        TOP SECRET                          %%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all;clc;
close all;
format compact;
format long;

Approximation_List = {'Butterworth','Chebyshev Type 1','Chebyshev Type 2'};
Approximation_Select = listdlg('ListString',Approximation_List','SelectionMode','single','ListSize',[250 75],'Name','Select Approximation:');

Filter_List = {'Low-Pass','High-Pass','Band-pass','Band-stop'};
Filter_Select = listdlg('ListString',Filter_List','SelectionMode','single','ListSize',[175 100],'Name','Select Filter:');

% ------------------------------------------------------------------------------------------------------------ %
Ap = str2double(cell2mat(inputdlg('Enter Pass-band Attenuation (Ap):')));
As = str2double(cell2mat(inputdlg('Enter Stop-band Attenuation (As):')));
if (Filter_Select==1 || Filter_Select==2)
    fp = str2double(cell2mat(inputdlg('Enter Pass-band cutoff Frequency (fp):')));
    fs = str2double(cell2mat(inputdlg('Enter Stop-band cutoff Frequency (fs):')));
end
if (Filter_Select==3)
    f1 = str2double(cell2mat(inputdlg('Enter Lower Stop-band cutoff Frequency (f1):')));
    f2 = str2double(cell2mat(inputdlg('Enter Lower Pass-band cutoff Frequency (f2):')));
    f3 = str2double(cell2mat(inputdlg('Enter Upper Pass-band cutoff Frequency (f3):')));
    f4 = str2double(cell2mat(inputdlg('Enter Upper Stop-band cutoff Frequency (f4):')));
end
if (Filter_Select==4)
    f1 = str2double(cell2mat(inputdlg('Enter Lower Pass-band cutoff Frequency (f1):')));
    f2 = str2double(cell2mat(inputdlg('Enter Lower Stop-band cutoff Frequency (f2):')));
    f3 = str2double(cell2mat(inputdlg('Enter Upper Stop-band cutoff Frequency (f3):')));
    f4 = str2double(cell2mat(inputdlg('Enter Upper Pass-band cutoff Frequency (f4):')));
end
fsampling = str2double(cell2mat(inputdlg('Enter sampling rate (fsampling):')));
excess = questdlg('Select Excess Tolerance','','Stopband','Passband','Stopband');
% ------------------------------------------------------------------------------------------------------------ %

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%% Buttersworth %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if (Approximation_Select == 1)
    %% Step 1 - Digital Frequencies ( f --> theta )
    if (Filter_Select==1 || Filter_Select==2)
        thetap = 2*pi*fp/fsampling;
        thetas = 2*pi*fs/fsampling;
    end
    if (Filter_Select==3 || Filter_Select==4)
        theta1 = 2*pi*f1/fsampling;
        theta2 = 2*pi*f2/fsampling;
        theta3 = 2*pi*f3/fsampling;
        theta4 = 2*pi*f4/fsampling;
    end
    
    %% Step 2 - Prewarp ( theta --> omega )
    if (Filter_Select==1 || Filter_Select==2)
        wp = 2*fsampling*tan(thetap/2);
        ws = 2*fsampling*tan(thetas/2);
    end
    if (Filter_Select==3 || Filter_Select==4)
        w1 = 2*fsampling*tan(theta1/2);
        w2 = 2*fsampling*tan(theta2/2);
        w3 = 2*fsampling*tan(theta3/2);
        w4 = 2*fsampling*tan(theta4/2);
    end
    
    %% Step 3 - Backware Transformation
    Wp = 1;
    if (Filter_Select==1)
        Ws = ws/wp;
    end
    if (Filter_Select==2)
        Ws = wp/ws;
    end
    if (Filter_Select==3)
        Ws1 = (w2*w3-w1^2)/(w1*(w3-w2));
        Ws2 = (w4^2-w2*w3)/(w4*(w3-w2));
        Ws = min(abs(Ws1),abs(Ws2));
        if (abs(Ws1) > abs(Ws2))
            w1a = (w2*w3)/w4;
            theta1a = 2*atan(w1a/(2*fsampling));
        else
            w4a = (w2*w3)/w1;
            theta4a = 2*atan(w4a/(2*fsampling));
        end
    end
    if (Filter_Select==4)
        Ws1 = (w2*(w4-w1))/(w1*w4-w2^2);
        Ws2 = (w3*(w4-w1))/((w3^2)-w1*w4);
        Ws = min(abs(Ws1),abs(Ws2));
        if (abs(Ws1) > abs(Ws2))
            w2a = (w1*w4)/w3;
            theta2a = 2*atan(w2a/(2*fsampling));
        else
            w3a = (w1*w4)/w2;
            theta3a = 2*atan(w3a/(2*fsampling));
        end
    end
    
    %% Step 4 - Order n of normalized LPF
    n = log10((10^(0.1*As)-1)/(10^(0.1*Ap)-1))/(2*log10(Ws));
    n = ceil(n);
    % Recalucation of Ap, As, epsilon, 3dB cutoff frequency Wc
    Ap1 = Ap;
    As1 = As;
    %if strcmp(excess,lower('Stopband'))                              %#ok<STCI>
    if excess == ('Stopband')                                         
        As1 = 10*log10(Ws^(2*n)*(10^(0.1*Ap1)-1)+1);
        epsilon = sqrt(10^(0.1*As1)-1)/Ws^n;
        Wc = 1/(epsilon^(1/n));
    else
        Ap1 = 10*log10((10^(0.1*As1)-1)/Ws^(2*n)+1);
        epsilon = sqrt(10^(0.1*Ap1)-1);
        Wc = 1/(epsilon^(1/n));
    end
    
    %% Step 5 - Normalized Lowpass Poles and Zeros
    k = 0:ceil((n-3)/2);
    S = exp(1i*pi*(2*k+n+1)/(2*n));
    S = [S,conj(S)];
    S = cplxpair(S);
    if rem(n,2)==1
        S=[-1,S];
    end
    if (Filter_Select==1 || Filter_Select==2)
        S = reshape(S,[1,n]);
    end
    SZ = realmax*ones(1,n);
    S = Wc*S;
    if (Filter_Select == 1 || Filter_Select == 2)
        wo = Wc*wp;
    end
    if (Filter_Select == 3 || Filter_Select == 4)
        wo = sqrt(w2*w3);
    end
    
    %% Step 6 - Frequency transformed poles and zeros
    if (Filter_Select==1)
        s = wp*S;
        sz = SZ;
    end
    if (Filter_Select==2)
        s = wp./S;
        sz = wp./SZ;
    end
    if (Filter_Select==3)
        s1 = (S.*(w3-w2)+sqrt(S.^2*((w3-w2)^2)-4*w2*w3))./2;
        s2 = (S.*(w3-w2)-sqrt(S.^2*((w3-w2)^2)-4*w2*w3))./2;
        s = [s1,s2];
        s = cplxpair(s);
        sz=[0*ones(1,n),realmax*ones(1,n)];
    end
    if (Filter_Select==4)
        s1 = ((w4-w1)+sqrt(((w4-w1)^2)-4*S.^2*w1*w4))./(2*S);
        s2 = ((w4-w1)-sqrt(((w4-w1)^2)-4*S.^2*w1*w4))./(2*S);
        s = [s1,s2];
        s = cplxpair(s);
        sz=1i*sqrt(w1*w4)*ones(1,n);
        sz=[sz,conj(sz)];
        sz=cplxpair(sz);
    end
    %% Step 7 - Bilinear Transformation
    pa = (2*fsampling+s)./(2*fsampling-s);
    za = (2*fsampling+sz)./(2*fsampling-sz);
    pa = cplxpair(pa);
    za = cplxpair(za);
    if rem(n,2)==1
        p = [pa,0];z = [za,0];
    else
        p = pa; z = za;
    end
    p;
    z;
    
    %% Step 8 - Second order sections and single section
    if (Filter_Select == 1 || Filter_Select == 2)
        p1 = reshape(p,[2,floor((n+1)/2)]);
        z1 = reshape(z,[2,floor((n+1)/2)]);
        DenSOS = [ones(1,floor((n+1)/2));-sum(p1);prod(p1)];
        NumSOS = [ones(1,floor((n+1)/2));-sum(z1);prod(z1)];
    end
    if (Filter_Select == 3)
        p1=reshape(pa,[2,n]);
        z1=reshape(za,[n,2]).';
        DenSOS = [ones(1,n);-sum(p1);prod(p1)];
        NumSOS = [ones(1,n);-sum(z1);prod(z1)];
    end
    if (Filter_Select == 4)
        p1=reshape(p,[2,n]);
        z1=reshape(z,[2,n]);
        DenSOS=[ones(1,n);-sum(p1);prod(p1)];
        NumSOS=[ones(1,n);-sum(z1);prod(z1)];
    end
    Den = poly(p);
    Num = poly(z);
    
    %% Step 9 - Constant Terms
    if (Filter_Select==1)
        b0 = sum(DenSOS)./sum(NumSOS);
        singleb0 = prod(b0);
    end
    if (Filter_Select==2)
        L = 0:(length(Den)-1);
        b0=abs([1 -1 1]*DenSOS)./abs([1 -1 1]*NumSOS);
        singleb0 = prod(b0);
    end
    if (Filter_Select==3)
        thetao = 2*atan(wo/(2*fsampling));
        Q = exp(-1i*(0:2)'*thetao);
        b0 = abs(DenSOS'*Q)./abs(NumSOS'*Q);
        singleb0 = prod(b0);
    end
    if (Filter_Select==4)
        b0=sum(DenSOS)./sum(NumSOS);
        singleb0=sum(Den)/sum(Num);
    end
    
    %% Plots
    x = -1:0.05:1;
    ty = sqrt(1-x.^2);
    figure(1);
    plot(x,ty,':b',x,-ty,':b');
    hold on;
    plot(real(pa),imag(pa),'bx','MarkerSize',12);hold on;
    plot(real(za),imag(za),'bo','MarkerSize',12);grid on;
    axis([-1,1,-1,1]);axis square;
    nstr = num2str(n);
    text(-0.9,0.1,nstr,'FontSize',12);
    title('Digital Pole-Zero Diagram');hold off;
    theta = 0:0.002:pi;
    Q = exp(-1i*(0:2)'*theta);
    H = singleb0*Num*exp(-1i*(0:length(Num)-1)'*theta)./(Den*exp(-1i*(0:length(Den)-1)'*theta));
    figure(2);plot(theta/pi,abs(H));grid on;
    xlabel('\theta/\pi');ylabel('|H(\theta)|');
    title('Magnitude Response (Linear)');
    Mag = 20*log10(abs(H));
    Th = -100;
    Mag = (Mag>=Th).*Mag+Th*(Mag<Th);
    figure(3);plot(theta/pi,Mag);grid on;
    xlabel('\theta/\pi');ylabel('20*log_1_0(|H(\theta)|)');
    title('Magnitude Response (dB)');
    figure(4);plot(theta/pi,180/pi*angle(H));grid on;
    xlabel('\theta/\pi');ylabel('\angleH(\theta) (degrees)');
    title('Phase Response');
    
    %% Step 10 - Verify
    H0 = abs(H(1));
    H1 = abs(H(length(H)));
    if (Filter_Select==1 || Filter_Select==2)
        thetao = 2*atan(wo/(2*fsampling));
        thetac = [thetap,thetas,thetao];
        Qc = exp(-1i*(0:2)'*thetac);
        Hc = singleb0*Num*exp(-1i*(0:length(Num)-1)'*thetac)./(Den*exp(-1i*(0:length(Den)-1)'*thetac));
        Hc = 20*log10(abs(Hc));
    end
    Ap1;
    As1;
    
    %% Reported Data
    
    format short;
    disp('Digital Cutoff Frequencies (rad) :');
    if (Filter_Select==1 || Filter_Select==2)
        thetap                                                          %#ok<*NOPTS>
        thetas
    end
    if (Filter_Select==3 || Filter_Select==4)
        theta1
        theta2
        theta3
        theta4
    end
    disp('-----------------------------------------------------------------------------');
    disp('Prewarped Cutoff Frequencies (rad/s) :');
    if (Filter_Select==1 || Filter_Select==2)
        wp
        ws
    end
    if (Filter_Select==3 || Filter_Select==4)
        w1
        w2
        w3
        w4
    end
    disp('-----------------------------------------------------------------------------');
    disp('Backward Cutoff Frequencies (rad/s) :');
    if (Filter_Select==1 || Filter_Select==2)
        Wp
        Ws
    end
    if (Filter_Select==3 || Filter_Select==4)
        Wp
        Ws1
        Ws2
        Ws
    end
    disp('-----------------------------------------------------------------------------');
    disp('Order ''n'' of the filter:');
    n
    disp('-----------------------------------------------------------------------------');
    disp('Normalized Poles:');
    S'
    disp('Normalized Zeros:');
    SZ'
    disp('-----------------------------------------------------------------------------');
    disp('Transformed Poles:');
    s'
    disp('Transformed Zeroes:');
    sz'
    disp('-----------------------------------------------------------------------------');
    disp('Digital Poles:');
    if (Filter_Select == 1 || Filter_Select==3)
        p'
    end
    if (Filter_Select == 2 || Filter_Select==4)
        pa'
    end
    disp('Digital Zeros:');
    if (Filter_Select == 1 || Filter_Select==3)
        z
    end
    if (Filter_Select == 2)
        za
    end
    if (Filter_Select == 4)
        za'
    end
    disp('-----------------------------------------------------------------------------');
    disp('Transfer Function Numerators of Digital Filter:');
    for x=1:length(NumSOS)
        NumSOS(:,x)'
    end
    disp('Transfer Function corresponding Denominators of Digital Filter:');
    for x=1:length(DenSOS)
        DenSOS(:,x)'
    end
    disp('First term is a constant, second term is z^-1, third term is z^-2');
    disp('-----------------------------------------------------------------------------');
    disp('b0 for each second-order section:');
    b0
    disp('Single b0:');
    singleb0
    disp('Verify (db):');
    Ap1
    As1
    disp('-----------------------------------------------------------------------------');

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%% Chebyshev Type 1 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if (Approximation_Select == 2)
    %% Step 1 - Digital Frequencies ( f --> theta )
    if (Filter_Select==1 || Filter_Select==2)
        thetap = 2*pi*fp/fsampling;
        thetas = 2*pi*fs/fsampling;
    end
    if (Filter_Select==3 || Filter_Select==4)
        theta1 = 2*pi*f1/fsampling;
        theta2 = 2*pi*f2/fsampling;
        theta3 = 2*pi*f3/fsampling;
        theta4 = 2*pi*f4/fsampling;
    end
    
    %% Step 2 - Prewarp ( theta --> omega )
    if (Filter_Select==1 || Filter_Select==2)
        wp = 2*fsampling*tan(thetap/2);
        ws = 2*fsampling*tan(thetas/2);
    end
    if (Filter_Select==3 || Filter_Select==4)
        w1 = 2*fsampling*tan(theta1/2);
        w2 = 2*fsampling*tan(theta2/2);
        w3 = 2*fsampling*tan(theta3/2);
        w4 = 2*fsampling*tan(theta4/2);
    end
    
    %% Step 3 - Backware Transformation
    Wp = 1;
    if (Filter_Select==1)
        Ws = ws/wp;
    end
    if (Filter_Select==2)
        Ws = wp/ws;
    end
    if (Filter_Select==3)
        Ws1 = (w2*w3-w1^2)/(w1*(w3-w2));
        Ws2 = (w4^2-w2*w3)/(w4*(w3-w2));
        Ws = min(abs(Ws1),abs(Ws2));
        if (abs(Ws1) > abs(Ws2))
            w1a = (w2*w3)/w4;
            theta1a = 2*atan(w1a/(2*fsampling));
        else
            w4a = (w2*w3)/w1;
            theta4a = 2*atan(w4a/(2*fsampling));
        end
    end
    if (Filter_Select==4)
        Ws1 = (w2*(w4-w1))/(w1*w4-w2^2);
        Ws2 = (w3*(w4-w1))/((w3^2)-w1*w4);
        Ws = min(abs(Ws1),abs(Ws2));
        if (abs(Ws1) > abs(Ws2))
            w2a = (w1*w4)/w3;
            theta2a = 2*atan(w2a/(2*fsampling));
        else
            w3a = (w1*w4)/w2;
            theta3a = 2*atan(w3a/(2*fsampling));
        end
    end
    
    %% Step 4 - Order n of normalized LPF
    n = acosh(sqrt((10^(0.1*As)-1)/(10^(0.1*Ap)-1)))/(acosh(Ws));
    n = ceil(n);
    % Recalucation of Ap, As, epsilon, 3dB cutoff frequency Wc
    Ap1 = Ap;
    As1 = As;
    %if strcmp(excess,lower('Stopband'))                              %#ok<STCI>
    if excess == ('Stopband')                                         
        As1 = 10*log10(1+(10^(0.1*Ap1)-1)*(cosh(n*acosh(Ws)))^2);
        epsilon = sqrt(10^(0.1*Ap1)-1);
    else
        Ap1 = 10*log10(1+(10^(0.1*As1)-1)/(cosh(n*acosh(Ws)))^2);
        epsilon = sqrt(10^(0.1*Ap1)-1);
    end
    
    %% Step 5 - Normalized Lowpass Poles and Zeros
    a = (1/n)*asinh(1/epsilon);
    k = 0:ceil((n-3)/2);
    S = -sinh(a)*sin((pi*(2*k+1))/(2*n)) + 1i*cosh(a)*cos((pi*(2*k+1))/(2*n));
    S = [S,conj(S)];
    S = cplxpair(S);
    if rem(n,2)==1
        S=[-sinh(a),S];
    end
    if (Filter_Select==1 || Filter_Select==2)
        S = reshape(S,[1,n]);
    end
    SZ = realmax*ones(1,n);
    
    %% Step 6 - Frequency transformed poles and zeros
    if (Filter_Select==1)
        s = wp*S;
        sz = SZ;
    end
    if (Filter_Select==2)
        s = wp./S;
        sz = wp./SZ;
    end
    if (Filter_Select==3)
        s1 = (S.*(w3-w2)+sqrt(S.^2*((w3-w2)^2)-4*w2*w3))./2;
        s2 = (S.*(w3-w2)-sqrt(S.^2*((w3-w2)^2)-4*w2*w3))./2;
        s = [s1,s2];
        s = cplxpair(s);
        sz=[0*ones(1,n),realmax*ones(1,n)];
    end
    if (Filter_Select==4)
        s1 = ((w4-w1)+sqrt(((w4-w1)^2)-4*S.^2*w1*w4))./(2*S);
        s2 = ((w4-w1)-sqrt(((w4-w1)^2)-4*S.^2*w1*w4))./(2*S);
        s = [s1,s2];
        s = cplxpair(s);
        sz=1i*sqrt(w1*w4)*ones(1,n);
        sz=[sz,conj(sz)];
        sz=cplxpair(sz);
    end
    %% Step 7 - Bilinear Transformation
    pa = (2*fsampling+s)./(2*fsampling-s);
    za = (2*fsampling+sz)./(2*fsampling-sz);
%     pa = ((2*fsampling./s)+1)./((2*fsampling./s)-1);
%     za = ((2*fsampling./sz)+1)./((2*fsampling./sz)-1);
    pa = cplxpair(pa);
    za = cplxpair(za);
    if (rem(n,2)==1 && Filter_Select ~= 4)
        za(1) = realmax*ones(1,1);
        p = [pa,0];z = [za,0];
    else
        p = pa; z = za;
    end
    p;
    z;
    
    %% Step 8 - Second order sections and single section
    if (Filter_Select == 1 || Filter_Select == 2)
        p1 = reshape(p,[2,floor((n+1)/2)]);
        z1 = reshape(z,[2,floor((n+1)/2)]);
        DenSOS = [ones(1,floor((n+1)/2));-sum(p1);prod(p1)];
        NumSOS = [ones(1,floor((n+1)/2));-sum(z1);prod(z1)];
    end
    if (Filter_Select == 3)
        p1=reshape(pa,[2,n]);
        z1=reshape(za,[n,2]).';
        DenSOS = [ones(1,n);-sum(p1);prod(p1)];
        NumSOS = [ones(1,n);-sum(z1);prod(z1)];
    end
    if (Filter_Select == 4)
        p1=reshape(p,[2,n]);
        z1=reshape(z,[2,n]);
        DenSOS=[ones(1,n);-sum(p1);prod(p1)];
        NumSOS=[ones(1,n);-sum(z1);prod(z1)];
    end
    Den = poly(p);
    Num = poly(z);
    
    %% Step 9 - Constant Terms
    if (Filter_Select==1)
        if (rem(n,2)==0)
            b0 = ((10^(-0.05*Ap1/(n/2)))*sum(DenSOS))./sum(NumSOS);
        else
            b0=sum(DenSOS)./sum(NumSOS);
        end
        singleb0 = prod(b0);
    end
    if (Filter_Select==2)
        if (rem(n,2)==0)
            b0 = abs(((10^(-0.05*Ap1/(n/2)))*([1 -1 1]*DenSOS)./abs([1 -1 1]*NumSOS)));
        else
            b0=abs([1 -1 1]*DenSOS)./abs([1 -1 1]*NumSOS);
        end
        singleb0 = prod(b0);       
    end
    if (Filter_Select==3)
        wo = sqrt(w2*w3);
        thetao = 2*atan(wo/(2*fsampling));
        Q = exp(-1i*(0:2)'*thetao);
        if (rem(n,2)==0)
            b0 = abs(((10^(-0.05*Ap1/(n/2)))*abs(DenSOS'*Q)./abs(NumSOS'*Q)));
        else
            b0 = abs(DenSOS'*Q)./abs(NumSOS'*Q);
        end
        singleb0 = prod(b0); 
    end
    if (Filter_Select==4)
        if (rem(n,2)==0)
            b0 = ((10^(-0.05*Ap1/(n/2)))*sum(DenSOS))./sum(NumSOS);
        else
            b0=sum(DenSOS)./sum(NumSOS);
        end
        singleb0 = prod(b0);
    end
    
    %% Plots
    x = -1:0.05:1;
    ty = sqrt(1-x.^2);
    figure(1);
    plot(x,ty,':b',x,-ty,':b');
    hold on;
    plot(real(pa),imag(pa),'bx','MarkerSize',12);hold on;
    plot(real(za),imag(za),'bo','MarkerSize',12);grid on;
    axis([-1,1,-1,1]);axis square;
    nstr = num2str(n);
    text(-0.9,0.1,nstr,'FontSize',12);
    title('Digital Pole-Zero Diagram');hold off;
    theta = 0:0.002:pi;
    Q = exp(-1i*(0:2)'*theta);
    H = singleb0*Num*exp(-1i*(0:length(Num)-1)'*theta)./(Den*exp(-1i*(0:length(Den)-1)'*theta));
    figure(2);plot(theta/pi,abs(H));grid on;
    xlabel('\theta/\pi');ylabel('|H(\theta)|');
    title('Magnitude Response (Linear)');
    Mag = 20*log10(abs(H));
    Th = -100;
    Mag = (Mag>=Th).*Mag+Th*(Mag<Th);
    figure(3);plot(theta/pi,Mag);grid on;
    xlabel('\theta/\pi');ylabel('20*log_1_0(|H(\theta)|)');
    title('Magnitude Response (dB)');
    figure(4);plot(theta/pi,180/pi*angle(H));grid on;
    xlabel('\theta/\pi');ylabel('\angleH(\theta) (degrees)');
    title('Phase Response');
    
    %% Step 10 - Verify
    H0 = abs(H(1));
    H1 = abs(H(length(H)));
    Ap1;
    As1;
    
    %% Reported Data
    
    format short;
    disp('Digital Cutoff Frequencies (rad) :');
    if (Filter_Select==1 || Filter_Select==2)
        thetap                                                          %#ok<*NOPTS>
        thetas
    end
    if (Filter_Select==3 || Filter_Select==4)
        theta1
        theta2
        theta3
        theta4
    end
    disp('-----------------------------------------------------------------------------');
    disp('Prewarped Cutoff Frequencies (rad/s) :');
    if (Filter_Select==1 || Filter_Select==2)
        wp
        ws
    end
    if (Filter_Select==3 || Filter_Select==4)
        w1
        w2
        w3
        w4
    end
    disp('-----------------------------------------------------------------------------');
    disp('Backward Cutoff Frequencies (rad/s) :');
    if (Filter_Select==1 || Filter_Select==2)
        Wp
        Ws
    end
    if (Filter_Select==3 || Filter_Select==4)
        Wp
        Ws1
        Ws2
        Ws
    end
    disp('-----------------------------------------------------------------------------');
    disp('Order ''n'' of the filter:');
    n
    disp('-----------------------------------------------------------------------------');
    disp('Normalized Poles:');
    S'
    disp('Normalized Zeros:');
    SZ'
    disp('-----------------------------------------------------------------------------');
    disp('Transformed Poles:');
    s'
    disp('Transformed Zeroes:');
    sz'
    disp('-----------------------------------------------------------------------------');
    disp('Digital Poles:');
    if (Filter_Select == 1 || Filter_Select==3)
        p'
    end
    if (Filter_Select == 2 || Filter_Select==4)
        pa'
    end
    disp('Digital Zeros:');
    if (Filter_Select == 1 || Filter_Select==3)
        z
    end
    if (Filter_Select == 2)
        za
    end
    if (Filter_Select == 4)
        za'
    end
    disp('-----------------------------------------------------------------------------');
    disp('Transfer Function Numerators of Digital Filter:');
    for x=1:length(NumSOS(1,:))
        NumSOS(:,x)'
    end
    disp('Transfer Function corresponding Denominators of Digital Filter:');
    for x=1:length(DenSOS(1,:))
        DenSOS(:,x)'
    end
    disp('First term is a constant, second term is z^-1, third term is z^-2');
    disp('-----------------------------------------------------------------------------');
    disp('b0 for each second-order section:');
    b0
    disp('Single b0:');
    singleb0
    disp('Verify (db):');
    Ap1
    As1
    disp('-----------------------------------------------------------------------------');

end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%% Chebyshev Type 2 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if (Approximation_Select == 3)
    %% Step 1 - Digital Frequencies ( f --> theta )
    if (Filter_Select==1 || Filter_Select==2)
        thetap = 2*pi*fp/fsampling;
        thetas = 2*pi*fs/fsampling;
    end
    if (Filter_Select==3 || Filter_Select==4)
        theta1 = 2*pi*f1/fsampling;
        theta2 = 2*pi*f2/fsampling;
        theta3 = 2*pi*f3/fsampling;
        theta4 = 2*pi*f4/fsampling;
    end
    
    %% Step 2 - Prewarp ( theta --> omega )
    if (Filter_Select==1 || Filter_Select==2)
        wp = 2*fsampling*tan(thetap/2);
        ws = 2*fsampling*tan(thetas/2);
    end
    if (Filter_Select==3 || Filter_Select==4)
        w1 = 2*fsampling*tan(theta1/2);
        w2 = 2*fsampling*tan(theta2/2);
        w3 = 2*fsampling*tan(theta3/2);
        w4 = 2*fsampling*tan(theta4/2);
    end
    
    %% Step 3 - Backware Transformation
    Wp = 1;
    if (Filter_Select==1)
        Ws = ws/wp;
    end
    if (Filter_Select==2)
        Ws = wp/ws;
    end
    if (Filter_Select==3)
        Ws1 = (w2*w3-w1^2)/(w1*(w3-w2));
        Ws2 = (w4^2-w2*w3)/(w4*(w3-w2));
        Ws = min(abs(Ws1),abs(Ws2));
        if (abs(Ws1) > abs(Ws2))
            w1a = (w2*w3)/w4;
            theta1a = 2*atan(w1a/(2*fsampling));
        else
            w4a = (w2*w3)/w1;
            theta4a = 2*atan(w4a/(2*fsampling));
        end
    end
    if (Filter_Select==4)
        Ws1 = (w2*(w4-w1))/(w1*w4-w2^2);
        Ws2 = (w3*(w4-w1))/((w3^2)-w1*w4);
        Ws = min(abs(Ws1),abs(Ws2));
        if (abs(Ws1) > abs(Ws2))
            w2a = (w1*w4)/w3;
            theta2a = 2*atan(w2a/(2*fsampling));
        else
            w3a = (w1*w4)/w2;
            theta3a = 2*atan(w3a/(2*fsampling));
        end
    end
    
    %% Step 4 - Order n of normalized LPF
    n = acosh(sqrt((10^(0.1*As)-1)/(10^(0.1*Ap)-1)))/(acosh(Ws));
    n = ceil(n);
    % Recalucation of Ap, As, epsilon, 3dB cutoff frequency Wc
    Ap1 = Ap;
    As1 = As;
    %if strcmp(excess,lower('Stopband'))                              %#ok<STCI>
    if excess == ('Stopband')                                         
        As1 = 10*log10(1+(10^(0.1*Ap1)-1)*(cosh(n*acosh(Ws)))^2);
        epsilon = 1/(sqrt(10^(0.1*As1)-1));
    else
        Ap1 = 10*log10(1+(10^(0.1*As1)-1)/(cosh(n*acosh(Ws)))^2);
        epsilon = 1/(sqrt(10^(0.1*As1)-1));
    end
    
    %% Step 5 - Normalized Lowpass Poles and Zeros
    a = (1/n)*asinh(sqrt(10^(0.1*As1)-1));
    k = 0:ceil((n-3)/2);
    S = -1./(sinh(a)*sin((pi*(2*k+1))/(2*n)) + 1i*cosh(a)*cos((pi*(2*k+1))/(2*n)));
    S = [S,conj(S)];
    S = cplxpair(S);
    if rem(n,2)==1
        S=[-1/sinh(a),S];
    end
    if (Filter_Select==1 || Filter_Select==2)
        S = reshape(S,[1,n]);
    end
    k1=0:floor((n-2)/2);
    SZ = 1i*sec((pi*(2*k1+1))/(2*n));
    SZ = [SZ,conj(SZ)];
    SZ = cplxpair(SZ);
    if rem(n,2)==1
        SZa=[inf,SZ];
    end
    if (Filter_Select==3)
        if excess == ('Stopband')                              %#ok<*BDSCA>
            Wc = cosh((1/n)*acosh(sqrt((10^(0.1*As1)-1)/(10^(0.1*Ap1)-1))));
        else
            Wc = (w4a-w1)/(w3-w2);
        end
        S = Wc*S;
        SZ = Wc*SZ;
    end
    if (Filter_Select==4)
        if excess == ('Stopband')                              %#ok<*BDSCA>
            Wc = cosh((1/n)*acosh(sqrt((10^(0.1*As1)-1)/(10^(0.1*Ap1)-1))));
        else
            Wc = (w4-w1)/(w3a-w2a);
        end
        S = Wc*S;
        SZ = Wc*SZ;
    end
    %% Step 6 - Frequency transformed poles and zeros
    if (Filter_Select==1)
        wo = wp*cosh((1/n)*acosh(sqrt((10^(0.1*As1)-1)/(10^(0.1*Ap1)-1))));
        if excess == ('Stopband')                              %#ok<*BDSCA>
            s = wo*S;
            sz = wo*SZ;
        else
            s = ws*S;
            sz = ws*SZ;
        end
    end
    if (Filter_Select==2)
        wo = wp/(cosh((1/n)*acosh(sqrt((10^(0.1*As1)-1)/(10^(0.1*Ap1)-1)))));
        if excess == ('Stopband')                              %#ok<*BDSCA>
            s = wo./S;
            sz = wo./SZ;
        else
            s = ws./S;
            sz = ws./SZ;
        end
    end
    if (Filter_Select==3)
        s1 = (S.*(w3-w2)+sqrt(S.^2*((w3-w2)^2)-4*w2*w3))./2;
        s2 = (S.*(w3-w2)-sqrt(S.^2*((w3-w2)^2)-4*w2*w3))./2;
        s = [s1,s2];
        s = cplxpair(s);
        sz1 = (SZ.*(w3-w2)+sqrt(SZ.^2*((w3-w2)^2)-4*w2*w3))./2;
        sz2 = (SZ.*(w3-w2)-sqrt(SZ.^2*((w3-w2)^2)-4*w2*w3))./2;
        sz = [sz1,sz2];
        sz = cplxpair(sz);
    end
    if (Filter_Select==4)
        s1 = ((w4-w1)+sqrt(((w4-w1)^2)-4*S.^2*w1*w4))./(2*S);
        s2 = ((w4-w1)-sqrt(((w4-w1)^2)-4*S.^2*w1*w4))./(2*S);
        s = [s1,s2];
        s = cplxpair(s);
        sz1 = 1i*sqrt(w1*w4);
        sz1 = [sz1,conj(sz1)];
        sz2 = ((w4-w1)+sqrt(((w4-w1)^2)-4*SZ.^2*w1*w4))./(2*SZ);
        sz3 = ((w4-w1)-sqrt(((w4-w1)^2)-4*SZ.^2*w1*w4))./(2*SZ);
        sz=[sz1,sz2,sz3];
        %sz=cplxpair(sz);
    end
    %% Step 7 - Bilinear Transformation
    pa = (2*fsampling+s)./(2*fsampling-s);
    za = (2*fsampling+sz)./(2*fsampling-sz);
    pa = cplxpair(pa);
    za = cplxpair(za);
%     if rem(n,2)==1
%         find(za==NaN)
%     end
%     za = cplxpair(za);
    if (rem(n,2)==1 && Filter_Select ~= 4)
        p = [pa,0];z = [za,-1];z2 = [za,-1,0];
    else
        p = pa; z2 = za; z = za;
    end
    
    %% Step 8 - Second order sections and single section
    if (Filter_Select == 1 || Filter_Select == 2)
        p1 = reshape(p,[2,floor((n+1)/2)]);
        z1 = reshape(z2,[2,floor((n+1)/2)]);
        DenSOS = [ones(1,floor((n+1)/2));-sum(p1);prod(p1)];
        NumSOS = [ones(1,floor((n+1)/2));-sum(z1);prod(z1)];
    end
    if (Filter_Select == 3)
        p1=reshape(p,[2,n]);
        z1=reshape(z2,[n,2]);
        DenSOS = [ones(1,n);-sum(p1);prod(p1)];
        NumSOS = [ones(1,n);-sum(z1);prod(z1)];
    end
    if (Filter_Select == 4)
        p1=reshape(p,[2,n]);
        z1=reshape(z,[2,n]);
        DenSOS=[ones(1,n);-sum(p1);prod(p1)];
        NumSOS=[ones(1,n);-sum(z1);prod(z1)];
    end
    Den = poly(p);
    Num = poly(z);
    
    %% Step 9 - Constant Terms
    if (Filter_Select==1)
        b0 = sum(DenSOS)./sum(NumSOS);
        singleb0 = prod(b0);
    end
    if (Filter_Select==2)
        L = 0:(length(Den)-1);
        b0=abs([1 -1 1]*DenSOS)./abs([1 -1 1]*NumSOS);
        singleb0 = prod(b0);
    end
    if (Filter_Select==3)
        wo = sqrt(w2*w3);
        thetao = 2*atan(wo/(2*fsampling));
        Q = exp(-1i*(0:2)'*thetao);
        b0 = abs(DenSOS'*Q)./abs(NumSOS'*Q);
        singleb0 = prod(b0); 
    end
    if (Filter_Select==4)
        b0=sum(DenSOS)./sum(NumSOS);
        singleb0 = prod(b0);
    end
    
    %% Plots
    x = -1:0.05:1;
    ty = sqrt(1-x.^2);
    figure(1);
    plot(x,ty,':b',x,-ty,':b');
    hold on;
    plot(real(pa),imag(pa),'bx','MarkerSize',12);hold on;
    plot(real(za),imag(za),'bo','MarkerSize',12);grid on;
    axis([-1,1,-1,1]);axis square;
    nstr = num2str(n);
    text(-0.9,0.1,nstr,'FontSize',12);
    title('Digital Pole-Zero Diagram');hold off;
    theta = 0:0.002:pi;
    Q = exp(-1i*(0:2)'*theta);
    H = singleb0*Num*exp(-1i*(0:length(Num)-1)'*theta)./(Den*exp(-1i*(0:length(Den)-1)'*theta));
    figure(2);plot(theta/pi,abs(H));grid on;
    xlabel('\theta/\pi');ylabel('|H(\theta)|');
    title('Magnitude Response (Linear)');
    Mag = 20*log10(abs(H));
    Th = -100;
    Mag = (Mag>=Th).*Mag+Th*(Mag<Th);
    figure(3);plot(theta/pi,Mag);grid on;
    xlabel('\theta/\pi');ylabel('20*log_1_0(|H(\theta)|)');
    title('Magnitude Response (dB)');
    figure(4);plot(theta/pi,180/pi*angle(H));grid on;
    xlabel('\theta/\pi');ylabel('\angleH(\theta) (degrees)');
    title('Phase Response');
    
    %% Step 10 - Verify
    H0 = abs(H(1));
    H1 = abs(H(length(H)));
    Ap1;
    As1;
    
    %% Reported Data
    
    format short;
    disp('Digital Cutoff Frequencies (rad) :');
    if (Filter_Select==1 || Filter_Select==2)
        thetap                                                          %#ok<*NOPTS>
        thetas
    end
    if (Filter_Select==3 || Filter_Select==4)
        theta1
        theta2
        theta3
        theta4
    end
    disp('-----------------------------------------------------------------------------');
    disp('Prewarped Cutoff Frequencies (rad/s) :');
    if (Filter_Select==1 || Filter_Select==2)
        wp
        ws
    end
    if (Filter_Select==3 || Filter_Select==4)
        w1
        w2
        w3
        w4
    end
    disp('-----------------------------------------------------------------------------');
    disp('Backward Cutoff Frequencies (rad/s) :');
    if (Filter_Select==1 || Filter_Select==2)
        Wp
        Ws
    end
    if (Filter_Select==3 || Filter_Select==4)
        Wp
        Ws1
        Ws2
        Ws
    end
    disp('-----------------------------------------------------------------------------');
    disp('Order ''n'' of the filter:');
    n
    disp('-----------------------------------------------------------------------------');
    disp('Normalized Poles:');
    S'
    disp('Normalized Zeros:');
%     if rem(n,2)==1
%         SZa'
%     else
        SZ'
% end
    disp('-----------------------------------------------------------------------------');
    disp('Transformed Poles:');
    s'
    disp('Transformed Zeroes:');
    if rem(n,2)==1
        sza = [inf,sz]'
    else
        sz'
    end
    disp('-----------------------------------------------------------------------------');
    disp('Digital Poles:');
    if (Filter_Select == 1 || Filter_Select==3)
        p'
    end
    if (Filter_Select == 2 || Filter_Select==4)
        pa'
    end
    disp('Digital Zeros:');
    if (Filter_Select == 1 || Filter_Select==3)
        z'
    end
    if (Filter_Select == 2)
        za'
    end
    if (Filter_Select == 4)
        za'
    end
    disp('-----------------------------------------------------------------------------');
    disp('Transfer Function Numerators of Digital Filter:');
    for x=1:length(NumSOS(1,:))
        NumSOS(:,x)'
    end
    disp('Transfer Function corresponding Denominators of Digital Filter:');
    for x=1:length(DenSOS(1,:))
        DenSOS(:,x)'
    end
    disp('First term is a constant, second term is z^-1, third term is z^-2');
    disp('-----------------------------------------------------------------------------');
    disp('b0 for each second-order section:');
    b0
    disp('Single b0:');
    singleb0
    disp('Verify (db):');
    Ap1
    As1
    disp('-----------------------------------------------------------------------------');

end

%     disp('Losses at cutoff frequencies (db):');
%     if (Filter_Select==1 || Filter_Select==2)
%         xaxis = theta/pi;
%         index_p = thetap/pi;
%         index_s = thetas/pi;
%         loc_p = min(find(xaxis>index_p));                               %#ok<*MXFND>
%         loc_s = min(find(xaxis>index_s));
%         mag_p = 20*log10(abs(H(loc_p)));
%         mag_s = 20*log10(abs(H(loc_s)));
%         Losses_p = abs(Ap-abs(mag_p))
%         Losses_s = abs(As-abs(mag_s))
%     end
%     if (Filter_Select==2 || Filter_Select==3)
%         xaxis = theta/pi;
%         index_1 = theta1/pi;
%         index_2 = theta2/pi;
%         index_3 = theta3/pi;
%         index_4 = theta4/pi;
%         loc_1 = min(find(xaxis>index_1));                               %#ok<*MXFND>
%         loc_2 = min(find(xaxis>index_2));
%         loc_3 = min(find(xaxis>index_3));
%         loc_4 = min(find(xaxis>index_4));
%         mag_1 = 20*log10(abs(H(loc_1)));
%         mag_2 = 20*log10(abs(H(loc_2)));
%         mag_3 = 20*log10(abs(H(loc_3)));
%         mag_4 = 20*log10(abs(H(loc_4)));
%         Losses_1 = abs(Ap-abs(mag_1))
%         Losses_2 = abs(As-abs(mag_2))
%         Losses_3 = abs(As-abs(mag_3))
%         Losses_4 = abs(As-abs(mag_4))
%     end