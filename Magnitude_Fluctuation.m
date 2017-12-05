% magnitude fluctuation computation

dt = 0.01;
for n = 10:10:30
    N = n^2;
    IR = ray_tracing(N,fs);
    E =  Energy_resp(IR, dt, fs);
    if temp == 0
        temp = E;
    end
    dev = dev + E-temp;
    temp = E;
end

dev = dev/3;