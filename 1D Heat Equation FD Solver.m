disp('Please choose 1 or 2: ')
disp('To Enter a Function for Intial data.........1');
disp('To Enter discrete points as Initial data....2');
option1 = input('Your choice: ');

%s = 0.25;
maxstep = 0;
msd = 0;

if (option1 == 1)
    disp('Enter the initial temperature function as y(x) at t=0 as string variable: ');
    disp('For example, y=''sin(x)''')
    disp('Then leave the input mode- type ''return'' and <enter>')
    keyboard
    disp('Enter: ');
    x1 = input('Xmin = ');
    x2 = input('Xmax = ');
    N = input('Total number of spatial nodes: ');
    t = input('Final time: ');
    M = input('Total number of time steps: ');
    tol = input('Tolerance: ');
    hs=input('Value of the constant heat source: ');
    h = (x2-x1)/(N-1);
    k = t/(M-1);
    s = k/(h*h);
    for j=1:M
        for i=1:N
            T(i,j)=0;
            if (j == 1)
                x = x1 + (i-1)*h;
                T(i,1) = eval(y);
            end
        end
    end
end
if (option1 == 2)
    u = input('Enter spatial nodes as row vector: ');
    v = input('Enter initial temperatures as row vector: ');
    display('Enter: ')
    t = input('Final time: ');
    M = input('Total number of time steps: ');
    tol = input('Tolerance: ');
    [Nc,N] = size(u);
    h = u(2) - u(1);
    k = t/(M-1);
    s = k/(h*h);
    for j=1:M
        for i=1:N
            T(i,j)=0;
            if (j == 1)
                T(i,1) = v(i);
            end
        end
    end
end

disp('Enter the heat source function hs(x) as string variables: )');
disp('For example hs= ''x^2''');
disp('Then leave the input mode, type <return> and enter');
keyboard


for i=1:N
    x=x1+ (i-1)*h;
    V(i,1)=k*eval(hs);
end

pause(1.50)

disp('Please choose 1, 2 or 3: ')
disp('Forward Difference in Time, Center in Space....1');
disp('Backward Difference in Time, Center in Space...2');
disp('Crank-Nicolson Method..........................3');
option2 = input('Your choice: ' );

if(option2 == 1)
    
   for i=1:N
       for j=1:N
           A(i,j)=0;
           if(i==j )
               A(i,j) = 1-2*s;
           end
           if(abs(i-j)==1)
               A(i,j) = s;
           end
       end
   end
   
%A(1,1) = 1;
%A(N,N) = 1;
%A(1,2) = 0;
%A(N,N-1) = 0;
  
   for j=1:M-1
       delmax = 0;
       T(:,j+1)=A*T(:,j)+ V(:,1);
       for i=1:N
           delta = abs(T(i,j+1)-T(i,j));
           if(delta > delmax)
               delmax = delta;
           end
       end
       if(delmax < tol)
           maxstep = j+1;
           break
       end
   end
end

if(option2 == 2)
    
     for i=1:N
       for j=1:N
           A(i,j)=0;
           if(i==j)
               A(i,j) = 1+2*s;
           end
           if(abs(i-j)==1)
               A(i,j) = (-1)*s;
           end
       end
     end
     
%A(1,1) = 1;
%A(N,N) = 1;
%A(1,2) = 0;
%A(N,N-1) = 0;
  
   for j=1:M-1
       delmax = 0;
       T(:,j+1) = A \ (T(:,j)+ V(:,1));
       for i=1:N
           delta = abs(T(i,j+1)-T(i,j));
           if(delta > delmax)
               delmax = delta;
           end
       end
       if(delmax < tol)
           maxstep = j+1;
           break
       end
   end
end

if(option2 == 3)
    
      for i=1:N
       for j=1:N
           A(i,j)=0;
           if(i==j)
               A(i,j) = 1+s;
           end
           if(abs(i-j)==1)
               A(i,j) = (-0.5)*s;
           end
       end
      end
   
      for i=1:N
       for j=1:N
           B(i,j)=0;
           if(i==j)
               B(i,j) = 1-s;
           end
           if(abs(i-j)==1)
               B(i,j) = 0.5*s;
           end
       end
      end
   
%A(1,1) = 1;
%A(N,N) = 1;
%A(1,2) = 0;
%A(N,N-1) = 0;

%B(1,1) = 1;
%B(N,N) = 1;
%B(1,2) = 0;
%B(N,N-1) = 0;
  
   for j=1:M-1
       delmax = 0;
       T(:,j+1) = (A \ B)* T(:,j) + A \ V(:,1);
       for i=1:N
           delta = abs(T(i,j+1)-T(i,j));
           if(delta > delmax)
               delmax = delta;
           end
       end
       if(delmax < tol)
           maxstep = j+1;
           break
       end
   end
end
maxstep = j+1;

for i=1:M
    t_axis(i) = k*(i-1);
    if(i<maxstep)
        rmsd_axis(i) = k*(i-1);
    end
end

for i=1:N
    x_axis(i) = x1 + h*(i-1);
end

for j=1:maxstep-1
    msd = 0;
    for i=1:N
       msd = msd + (T(i,j+1)-T(i,j))^2;
    end
    msd = msd/N;
    rmsd(j) = sqrt(msd);
end
   
disp('Please choose: ');
disp('For simulation of Temperature against x-axis................1');
disp('For surface curve of Temperature against time and x-axis....2');
disp('For RMSD curve against time.................................3');

option3 = input('Your choice: ');

if( option3==1 )

for i=1:maxstep
    plot(x_axis,T(:,i),'r')
    ylim([min(T(:,1)) max(T(:,1))])
    grid on;
    xlabel('X-axis ->');
    ylabel('Temperature ->');
    title('Temperature Vs X-axis');
    pause(0.005);
    %hold on;
end
end

if( option3==2)
    
    surf(t_axis,x_axis,T)

end

if( option3==3)
    
    plot(rmsd_axis,rmsd,'r')
    grid on;
    xlabel('Time ->');
    ylabel('RMS Difference Average ->');
    title('');
    
end

fprintf('Final temperature profile: %d \n',T(:,maxstep));
plot(x_axis,T(:,maxstep))
fprintf('Maximum difference at final step: %d \n',delmax);
fprintf('RMS difference at final step: %d \n',rmsd(maxstep-1));
fprintf('Total times steps taken: %d \n',maxstep);


if (maxstep==M && delmax>tol)
    disp('To achieve given tolerance, increase number of time steps!');
end

