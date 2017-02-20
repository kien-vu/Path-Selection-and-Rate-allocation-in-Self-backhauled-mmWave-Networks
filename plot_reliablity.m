


network_queue0 = network_queue./(alpha1/4);
total_data = nnz(network_queue0(:,:,2:Iterations)>0);
up = max(max(max(network_queue0)))
N0 = round(up)/0.02;
A = zeros(N0,1);
a = zeros(N0,1);
a(1,1) = 0.02;
for n = 1:N0
    A(n,1) = (total_data - nnz(network_queue0(:,:,2:Iterations) >= a(n,1)))/total_data;
    a(n+1,1) = a(n,1) + 0.02;
end


network_queue1 = network_queue./(alpha1/2);
total_data = nnz(network_queue1(:,:,2:Iterations)>0);
up1 = max(max(max(network_queue1)));
N1 = round(up1)/0.02;
B = zeros(N1,1);
b = zeros(N1,1);
b(1,1) = 0.02;
for n = 1:N1
    B(n,1) = (total_data - nnz(network_queue1(:,:,2:Iterations) >= b(n,1)))/total_data;
    b(n+1,1) = b(n,1) + 0.02;
end

network_queue2 = network_queue./(alpha1/2);
total_data = nnz(network_queue2(:,:,2:Iterations)>0);
up2 = max(max(max(network_queue2)));
N2 = round(up2)/0.02;
C = zeros(N2,1);
c = zeros(N2,1);
 c(1,1) = 0.02;
for n = 1:N2
    C(n,1) = (total_data - nnz(network_queue2(:,:,2:Iterations) >= c(n,1)))/total_data;
    c(n+1,1) = c(n,1) + 0.02;
end

figure
plot(2.5*c(1:N2), C); hold on


figure
% plot(a(1:N0), A); hold on
plot(c(1:N2), C); hold on
plot(b(1:N1), B); hold on
xlabel('Latency requirement d^{th} [ms]','FontSize',12.6,'FontName','Times New Roman')
ylabel('Probability that Pr(delay \leq d^{th}) = 1- \epsilon ','FontSize',12.6,'FontName','Times New Roman')

for n = 1:N2
% while(1)
        if C(n)>=0.95%(5.9 <= c(n)) && (c(n) <= 6.1)
            break;
        end
% end
end