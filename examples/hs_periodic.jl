using LinearAlgebra
using Plots
using RecursiveArrayTools
using Revise
using Uchiyama

function main( nstep )
    N = 100
    ϵ = 0.02
    
    q = VectorOfArray([zeros(2) for i in 1:N])
    q[1] = [0.5,0.5]
    
    for klm = 2:N
        frein = 0
        overlap = 1
        p = zeros(2)
        while ((overlap == 1) && frein < 7000)
            p .= 2ϵ * ones(2) + (1 - 4ϵ) * rand(2)
            overlap2 = 0
            for subh = 1:(klm-1)
                if norm(p - q[subh],2) < 2ϵ
                    overlap2 = 1
                end
            end
            overlap = overlap2
            frein += 1
            if frein == 10000
                frein
                @error "frein trop grand"
                exit
            end
        end
        q[klm] = p
    end
    
    
    offset = [[0, 0], [1,  0], [-1,  0], 
              [0, 1], [0, -1], [-1,  1], 
              [1, 1], [1, -1], [-1, -1]]
    
    function compute_dt( qi, qj, vi, vj )
    
        δz = qi .- qj
        δv = vi .- vj
    
        c = (δv'δz).^2 - (δv'δv) * (δz'δz - (2ϵ)^2)
    
        if δv'δz >= 0
            δt = Inf
        elseif c < 0
            δt = Inf
        else
            δt = -(δv'δz + sqrt(c)) / (δv'δv)
        end
    
        return δt
    
    end
    
    
    v = VectorOfArray([randn(2) for j = 1:N])
    
    Collisions = Inf .* ones(N,N)
    Fantome = zeros(Int, N,N)
    
    for k=1:N
        for l=k+1:N
            t1 = compute_dt(q[l]+[ 0, 0], q[k], v[l], v[k])
            t2 = compute_dt(q[l]+[ 1, 0], q[k], v[l], v[k])
            t3 = compute_dt(q[l]+[-1, 0], q[k], v[l], v[k])
            t4 = compute_dt(q[l]+[ 0, 1], q[k], v[l], v[k])
            t5 = compute_dt(q[l]+[ 0,-1], q[k], v[l], v[k])
            t6 = compute_dt(q[l]+[-1, 1], q[k], v[l], v[k]) 
            t7 = compute_dt(q[l]+[ 1, 1], q[k], v[l], v[k])
            t8 = compute_dt(q[l]+[ 1,-1], q[k], v[l], v[k])
            t9 = compute_dt(q[l]+[-1,-1], q[k], v[l], v[k])
            a  = [t1,t2,t3,t4,t5,t6,t7,t8,t9]
            Collisions[k,l] = minimum(a)
            if isinf(Collisions[k,l])
                Fantome[k,l] = 0
            else
                Fantome[k,l] = argmin(a)
            end
        end
    end
    
    anim = @animate for step in 1:nstep
    
        m, n = Tuple(argmin(Collisions))
        dt = Collisions[m,n]
        
        if Fantome[m,n] == 1  
    
            qa = [q[n] + offset[k] + dt*v[n] for k in 1:9]
            qb = q[m] + dt .* v[m]
            vr2 = argmin([norm(qx-qb) for qx in qa])
    
            q[n]=qa[vr2]
            q[m]=qb
    
            J = (dot(v[n]-v[m],q[n]-q[m]))/(2ϵ)
            v[1,m]=v[1,m] + J* (q[1,n]-q[1,m])/(2ϵ)
            v[2,m]=v[2,m] + J* (q[2,n]-q[2,m])/(2ϵ)
            v[1,n]=v[1,n] - J* (q[1,n]-q[1,m])/(2ϵ)
            v[2,n]=v[2,n] - J* (q[2,n]-q[2,m])/(2ϵ)
            q[n] = mod.(q[n],1)
            q[m] = mod.(q[m],1)
         
         elseif Fantome[m,n] == 2 
    
            qa = [q[n] + offset[k] + dt*v[n] for k in 1:9]
             qa1 = q[n] + [ 0, 0] + dt*v[n]
             qa2 = q[n] + [ 1, 0] + dt*v[n]
             qa3 = q[n] + [-1, 0] + dt*v[n]
             qa4 = q[n] + [ 0, 1] + dt*v[n]
             qa5 = q[n] + [ 0,-1] + dt*v[n]
             qa6 = q[n] + [-1, 1] + dt*v[n]
             qa7 = q[n] + [ 1, 1] + dt*v[n]
             qa8 = q[n] + [ 1,-1] + dt*v[n]
             qa9 = q[n] + [-1;-1] + dt*v[n]
             qb  = q[m] + dt*v[m]
    
             vraicoll=[norm(qa1-qb),norm(qa2-qb),norm(qa3-qb),
                       norm(qa4-qb),norm(qa5-qb),norm(qa6-qb),
                       norm(qa7-qb),norm(qa8-qb),norm(qa9-qb)]
    
             vr2 = argmin(vraicoll)
    
             if vr2==1
                 q[n]=qa1   
                 q[m]=qb
             elseif vr2 ==2
                 q[n]=qa2   
                 q[m]=qb
             elseif vr2==3
                 q[n]=qa3   
                 q[m]=qb
             elseif vr2==4
                 q[n]=qa4   
                 q[m]=qb
             elseif vr2==5
                 q[n]=qa5   
                 q[m]=qb
             elseif vr2==6
                 q[n]=qa6   
                 q[m]=qb
             elseif vr2==7
                 q[n]=qa7   
                 q[m]=qb
             elseif vr2==8
                 q[n]=qa8   
                 q[m]=qb
             elseif vr2==9
                 q[n]=qa9   
                 q[m]=qb 
             end
             J= (dot(v[n]-v[m],q[n]-q[m]))/(2ϵ)
             v[1,m]=v[1,m] + J* (q[1,n] - q[1,m])/(2ϵ)
             v[2,m]=v[2,m] + J* (q[2,n] - q[2,m])/(2ϵ)
             v[1,n]=v[1,n] - J* (q[1,n] - q[1,m])/(2ϵ)
             v[2,n]=v[2,n] - J* (q[2,n] - q[2,m])/(2ϵ)
             q[n]=mod.(q[n],1)
             q[m]=mod.(q[m],1)
         
         elseif Fantome[m,n] == 3 
     
            qa = [q[n] + offset[k] + dt*v[n] for k in 1:9]
             qa1= q[n] + [ 0, 0] + dt*v[n]
             qa2= q[n] + [ 1, 0] + dt*v[n]
             qa3= q[n] + [-1, 0] + dt*v[n]
             qa4= q[n] + [ 0, 1] + dt*v[n]
             qa5= q[n] + [ 0,-1] + dt*v[n]
             qa6= q[n] + [-1, 1] + dt*v[n]
             qa7= q[n] + [ 1, 1] + dt*v[n]
             qa8= q[n] + [ 1,-1] + dt*v[n]
             qa9= q[n] + [-1,-1] + dt*v[n]
             qb = q[m] + dt*v[m]
    
             vraicoll=[norm(qa1-qb), norm(qa2-qb), norm(qa3-qb),
                       norm(qa4-qb), norm(qa5-qb), norm(qa6-qb),
                       norm(qa7-qb), norm(qa8-qb), norm(qa9-qb)]
    
             vr2 = argmin(vraicoll)
    
             if vr2==1
                 q[n]=qa1   
                 q[m]=qb
             elseif vr2 ==2
                 q[n]=qa2   
                 q[m]=qb
             elseif vr2==3
                 q[n]=qa3   
                 q[m]=qb
             elseif vr2==4
                 q[n]=qa4   
                 q[m]=qb
             elseif vr2==5
                 q[n]=qa5   
                 q[m]=qb
              elseif vr2==6
                 q[n]=qa6   
                 q[m]=qb
              elseif vr2==7
                 q[n]=qa7   
                 q[m]=qb
              elseif vr2==8
                 q[n]=qa8   
                 q[m]=qb
              elseif vr2==9
                 q[n]=qa9   
                 q[m]=qb 
             end
    
             J= (dot(v[n]-v[m],q[n]  -q[m]))/(2ϵ)
    
             v[1,m]=v[1,m] + J* (q[1,n] -q[1,m])/(2ϵ)
             v[2,m]=v[2,m] + J* (q[2,n] -q[2,m])/(2ϵ)
             v[1,n]=v[1,n] - J* (q[1,n]  -q[1,m])/(2ϵ)
             v[2,n]=v[2,n] - J* (q[2,n]  -q[2,m])/(2ϵ)
    
             q[n]=mod.(q[n],1)
             q[m]=mod.(q[m],1)
         
         elseif Fantome[m,n] == 4
     
            qa = [q[n] + offset[k] + dt*v[n] for k in 1:9]
             qa1= q[n] + dt*v[n]
             qa2= q[n] + [ 1, 0] + dt*v[n]
             qa3= q[n] + [-1, 0] + dt*v[n]
             qa4= q[n] + [ 0, 1] + dt*v[n]
             qa5= q[n] + [ 0,-1] + dt*v[n]
             qa6= q[n] + [-1, 1] + dt*v[n]
             qa7= q[n] + [ 1, 1] + dt*v[n]
             qa8= q[n] + [ 1,-1] + dt*v[n]
             qa9= q[n] + [-1,-1] + dt*v[n]
             qb= q[m]  + dt*v[m]
             vraicoll=[norm(qa1-qb),norm(qa2-qb),norm(qa3-qb),
                       norm(qa4-qb),norm(qa5-qb),norm(qa6-qb),
                       norm(qa7-qb),norm(qa8-qb),norm(qa9-qb)]
            
             vr2 = argmin(vraicoll)
    
             if vr2==1
                 q[n]=qa1   
                 q[m]=qb
             elseif vr2 ==2
                 q[n]=qa2   
                 q[m]=qb
             elseif vr2==3
                 q[n]=qa3   
                 q[m]=qb
             elseif vr2==4
                 q[n]=qa4   
                 q[m]=qb
             elseif vr2==5
                 q[n]=qa5   
                 q[m]=qb
             elseif vr2==6
                 q[n]=qa6   
                 q[m]=qb
             elseif vr2==7
                 q[n]=qa7   
                 q[m]=qb
             elseif vr2==8
                 q[n]=qa8   
                 q[m]=qb
             elseif vr2==9
                 q[n]=qa9   
                 q[m]=qb 
             end
    
             J= (dot(v[n]-v[m],q[n]  -q[m]))/(2ϵ)
             v[1,m]=v[1,m] + J* (q[1,n] -q[1,m])/(2ϵ)
             v[2,m]=v[2,m] + J* (q[2,n] -q[2,m])/(2ϵ)
             v[1,n]=v[1,n] - J* (q[1,n] -q[1,m])/(2ϵ)
             v[2,n]=v[2,n] - J* (q[2,n] -q[2,m])/(2ϵ)
             q[n]=mod.(q[n],1)
             q[m]=mod.(q[m],1)
         
         elseif Fantome[m,n] == 5
     
            qa = [q[n] + offset[k] + dt*v[n] for k in 1:9]
             qa1= q[n]+ [ 0, 0] + dt*v[n]
             qa2= q[n]+ [ 1, 0] + dt*v[n]
             qa3= q[n]+ [-1, 0] + dt*v[n]
             qa4= q[n]+ [ 0, 1] + dt*v[n]
             qa5= q[n]+ [ 0,-1] + dt*v[n]
             qa6= q[n]+ [-1, 1] + dt*v[n]
             qa7= q[n]+ [ 1, 1] + dt*v[n]
             qa8= q[n]+ [ 1,-1] + dt*v[n]
             qa9= q[n]+ [-1,-1] + dt*v[n]
             qb = q[m] + dt*v[m]
             vraicoll=[norm(qa1-qb,2),norm(qa2-qb,2),norm(qa3-qb,2),
                       norm(qa4-qb,2),norm(qa5-qb,2),norm(qa6-qb,2),
                       norm(qa7-qb,2),norm(qa8-qb,2),norm(qa9-qb,2)]
    
             vr2 = argmin(vraicoll)
    
                if vr2==1
                 q[n]=qa1   
                 q[m]=qb
             elseif vr2 ==2
                 q[n]=qa2   
                 q[m]=qb
             elseif vr2==3
                 q[n]=qa3   
                 q[m]=qb
             elseif vr2==4
                 q[n]=qa4   
                 q[m]=qb
             elseif vr2==5
                 q[n]=qa5   
                 q[m]=qb
              elseif vr2==6
                 q[n]=qa6   
                 q[m]=qb
              elseif vr2==7
                 q[n]=qa7   
                 q[m]=qb
              elseif vr2==8
                 q[n]=qa8   
                 q[m]=qb
              elseif vr2==9
                 q[n]=qa9   
                 q[m]=qb 
             end
             J= (dot(v[n]-v[m],q[n]  -q[m]))/(2ϵ)
             v[1,m] = v[1,m] + J* (q[1,n] -q[1,m])/(2ϵ)
             v[2,m] = v[2,m] + J* (q[2,n] -q[2,m])/(2ϵ)
             v[1,n] = v[1,n] - J* (q[1,n] -q[1,m])/(2ϵ)
             v[2,n] = v[2,n] - J* (q[2,n] -q[2,m])/(2ϵ)
             q[n]=mod.(q[n],1)
             q[m]=mod.(q[m],1)
         
          elseif Fantome[m,n] == 6
     
            qa = [q[n] + offset[k] + dt*v[n] for k in 1:9]
              qa1= q[n]+ [ 0, 0] + dt*v[n]
              qa2= q[n]+ [ 1, 0] + dt*v[n]
              qa3= q[n]+ [-1, 0] + dt*v[n]
              qa4= q[n]+ [ 0, 1] + dt*v[n]
              qa5= q[n]+ [ 0,-1] + dt*v[n]
              qa6= q[n]+ [-1, 1] + dt*v[n]
              qa7= q[n]+ [ 1, 1] + dt*v[n]
              qa8= q[n]+ [ 1,-1] + dt*v[n]
              qa9= q[n]+ [-1,-1] + dt*v[n]
    
              qb= q[m] + dt*v[m]
    
              vraicoll=[norm(qa1-qb),norm(qa2-qb),norm(qa3-qb),
                        norm(qa4-qb),norm(qa5-qb),norm(qa6-qb),
                        norm(qa7-qb),norm(qa8-qb),norm(qa9-qb)]
    
              vr2 = argmin(vraicoll)
    
                  if vr2==1
                  q[n]=qa1   
                  q[m]=qb
              elseif vr2 ==2
                  q[n]=qa2   
                  q[m]=qb
              elseif vr2==3
                  q[n]=qa3   
                  q[m]=qb
              elseif vr2==4
                  q[n]=qa4   
                  q[m]=qb
              elseif vr2==5
                  q[n]=qa5   
                  q[m]=qb
              elseif vr2==6
                  q[n]=qa6   
                  q[m]=qb
              elseif vr2==7
                  q[n]=qa7   
                  q[m]=qb
              elseif vr2==8
                  q[n]=qa8   
                  q[m]=qb
              elseif vr2==9
                  q[n]=qa9   
                  q[m]=qb 
              end
              J= (dot(v[n]-v[m],q[n]  -q[m]))/(2ϵ)
              v[1,m]=v[1,m] + J* (q[1,n] -q[1,m])/(2ϵ)
              v[2,m]=v[2,m] + J* (q[2,n] -q[2,m])/(2ϵ)
              v[1,n]=v[1,n] - J* (q[1,n] -q[1,m])/(2ϵ)
              v[2,n]=v[2,n] - J* (q[2,n] -q[2,m])/(2ϵ)
              q[n]=mod.(q[n],1)
              q[m]=mod.(q[m],1)
         
         elseif Fantome[m,n] == 7 
     
            qa = [q[n] + offset[k] + dt*v[n] for k in 1:9]
             qa1= q[n]+ [ 0, 0] + dt*v[n]
             qa2= q[n]+ [ 1, 0] + dt*v[n]
             qa3= q[n]+ [-1, 0] + dt*v[n]
             qa4= q[n]+ [ 0, 1] + dt*v[n]
             qa5= q[n]+ [ 0,-1] + dt*v[n]
             qa6= q[n]+ [-1, 1] + dt*v[n]
             qa7= q[n]+ [ 1, 1] + dt*v[n]
             qa8= q[n]+ [ 1,-1] + dt*v[n]
             qa9= q[n]+ [-1,-1] + dt*v[n]
             qb= q[m] + dt*v[m]
             vraicoll=[norm(qa1-qb,2),norm(qa2-qb,2),norm(qa3-qb,2),
                       norm(qa4-qb,2),norm(qa5-qb,2),norm(qa6-qb,2),
                       norm(qa7-qb,2),norm(qa8-qb,2),norm(qa9-qb,2)]
            
             vr2 = argmin(vraicoll)
    
             if vr2==1
                 q[n]=qa1   
                 q[m]=qb
             elseif vr2 ==2
                 q[n]=qa2   
                 q[m]=qb
             elseif vr2==3
                 q[n]=qa3   
                 q[m]=qb
             elseif vr2==4
                 q[n]=qa4   
                 q[m]=qb
             elseif vr2==5
                 q[n]=qa5   
                 q[m]=qb
             elseif vr2==6
                 q[n]=qa6   
                 q[m]=qb
             elseif vr2==7
                 q[n]=qa7   
                 q[m]=qb
             elseif vr2==8
                 q[n]=qa8   
                 q[m]=qb
             elseif vr2==9
                 q[n]=qa9   
                 q[m]=qb 
             end
             J= (dot(v[n]-v[m],q[n]-q[m]))/(2ϵ)
             v[1,m]=v[1,m] + J* (q[1,n]  -q[1,m])/(2ϵ)
             v[2,m]=v[2,m] + J* (q[2,n]   -q[2,m])/(2ϵ)
             v[1,n]=v[1,n] - J* (q[1,n] -q[1,m])/(2ϵ)
             v[2,n]=v[2,n] - J* (q[2,n]  -q[2,m])/(2ϵ)
             q[n]=mod.(q[n],1)
             q[m]=mod.(q[m],1)
         
         elseif Fantome[m,n] == 8    
     
            qa = [q[n] + offset[k] + dt*v[n] for k in 1:9]
             qa1= q[n]+ [ 0, 0] + dt*v[n]
             qa2= q[n]+ [ 1, 0] + dt*v[n]
             qa3= q[n]+ [-1, 0] + dt*v[n]
             qa4= q[n]+ [ 0, 1] + dt*v[n]
             qa5= q[n]+ [ 0,-1] + dt*v[n]
             qa6= q[n]+ [-1, 1] + dt*v[n]
             qa7= q[n]+ [ 1, 1] + dt*v[n]
             qa8= q[n]+ [ 1,-1] + dt*v[n]
             qa9= q[n]+ [-1,-1] + dt*v[n]
             qb= q[m] + dt .* v[m]
             vraicoll=[norm(qa1-qb,2),norm(qa2-qb,2),norm(qa3-qb,2),
                       norm(qa4-qb,2),norm(qa5-qb,2),norm(qa6-qb,2),
                       norm(qa7-qb,2),norm(qa8-qb,2),norm(qa9-qb,2)]
    
             vr2 = argmin(vraicoll)
    
             if vr2==1
                 q[n]=qa1   
                 q[m]=qb
             elseif vr2 ==2
                 q[n]=qa2   
                 q[m]=qb
             elseif vr2==3
                 q[n]=qa3   
                 q[m]=qb
             elseif vr2==4
                 q[n]=qa4   
                 q[m]=qb
             elseif vr2==5
                 q[n]=qa5   
                 q[m]=qb
             elseif vr2==6
                 q[n]=qa6   
                 q[m]=qb
             elseif vr2==7
                 q[n]=qa7   
                 q[m]=qb
             elseif vr2==8
                 q[n]=qa8   
                 q[m]=qb
             elseif vr2==9
                 q[n]=qa9   
                 q[m]=qb 
             end
             J= (dot(v[n]-v[m],q[n]  -q[m]))/(2ϵ)
             v[1,m]=v[1,m] + J* (q[1,n]  -q[1,m])/(2ϵ)
             v[2,m]=v[2,m] + J* (q[2,n]  -q[2,m])/(2ϵ)
             v[1,n]=v[1,n] - J* (q[1,n]  -q[1,m])/(2ϵ)
             v[2,n]=v[2,n] - J* (q[2,n]  -q[2,m])/(2ϵ)
             q[n]=mod.(q[n],1)
             q[m]=mod.(q[m],1)
         
         elseif Fantome[m,n] == 9
     
            qa = [q[n] + offset[k] + dt*v[n] for k in 1:9]
             qa1= q[n]+ [ 0, 0] + dt*v[n]
             qa2= q[n]+ [ 1, 0] + dt*v[n]
             qa3= q[n]+ [-1, 0] + dt*v[n]
             qa4= q[n]+ [ 0, 1] + dt*v[n]
             qa5= q[n]+ [ 0,-1] + dt*v[n]
             qa6= q[n]+ [-1, 1] + dt*v[n]
             qa7= q[n]+ [ 1, 1] + dt*v[n]
             qa8= q[n]+ [ 1,-1] + dt*v[n]
             qa9= q[n]+ [-1,-1] + dt*v[n]
             qb= q[m] + dt*v[m]
             vraicoll=[norm(qa1-qb),norm(qa2-qb),norm(qa3-qb),
                       norm(qa4-qb),norm(qa5-qb),norm(qa6-qb),
                       norm(qa7-qb),norm(qa8-qb),norm(qa9-qb)]
    
             vr2 = argmin(vraicoll)
    
             if vr2==1
                 q[n]=qa1   
                 q[m]=qb
             elseif vr2 ==2
                 q[n]=qa2   
                 q[m]=qb
             elseif vr2==3
                 q[n]=qa3   
                 q[m]=qb
             elseif vr2==4
                 q[n]=qa4   
                 q[m]=qb
             elseif vr2==5
                 q[n]=qa5   
                 q[m]=qb
             elseif vr2==6
                 q[n]=qa6   
                 q[m]=qb
             elseif vr2==7
                 q[n]=qa7   
                 q[m]=qb
             elseif vr2==8
                 q[n]=qa8   
                 q[m]=qb
             elseif vr2==9
                 q[n]=qa9   
                 q[m]=qb 
             end
             J= (dot(v[n]-v[m],q[n]  -q[m]))/(2ϵ)
             v[1,m]=v[1,m] + J* (q[1,n] -q[1,m])/(2ϵ)
             v[2,m]=v[2,m] + J* (q[2,n] -q[2,m])/(2ϵ)
             v[1,n]=v[1,n] - J* (q[1,n] -q[1,m])/(2ϵ)
             v[2,n]=v[2,n] - J* (q[2,n] -q[2,m])/(2ϵ)
             q[n]=mod.(q[n],1)
             q[m]=mod.(q[m],1)
         
         end
         
         if m<n
            for i=1:m-1
                q[i]=q[i]+ dt*v[i]
                q[i]=mod.(q[i],1)
            end
            for i=m+1:n-1
                q[i]=q[i]+ dt*v[i]
                q[i]=mod.(q[i],1)
            end
            for i=n+1:N
                q[i]=q[i]+ dt*v[i]
                q[i]=mod.(q[i],1)
            end
        elseif n<m
            for i=1:n-1
                q[i]=q[i]+ dt*v[i]
                q[i]=mod.(q[i],1)
            end
            for i=n+1:m-1
                q[i]=q[i]+ dt*v[i]
                q[i]=mod.(q[i],1)
            end
            for i=m+1:N
                q[i]=q[i]+ dt*v[i]
                q[i]=mod.(q[i],1)
            end
        end

        plt = plot(size  = (200,200), 
                   xlims = (0,1), 
                   ylims = (0,1), 
                   grid  = false, 
                   axis  = nothing, 
                   framestyle= :none, 
                   legend=false, 
                   widen = false)
    
        c = trunc(Int, 200ϵ)
        scatter!( q[1,:], q[2,:], markershape  = :circle, 
              markersize   = c , 
              aspect_ratio = :equal)

         
        Collisions .-= dt
        Collisions[m,:] .= Inf
        Collisions[n,:] .= Inf
        Collisions[:,m] .= Inf
        Collisions[:,n] .= Inf
        
        for l=1:m-1
            t1=compute_dt(q[l]+[ 0, 0],q[m],v[l],v[m])
            t2=compute_dt(q[l]+[ 1, 0],q[m],v[l],v[m])
            t3=compute_dt(q[l]+[-1, 0],q[m],v[l],v[m])
            t4=compute_dt(q[l]+[ 0, 1],q[m],v[l],v[m])
            t5=compute_dt(q[l]+[ 0,-1],q[m],v[l],v[m])
            t6=compute_dt(q[l]+[-1, 1],q[m],v[l],v[m]) 
            t7=compute_dt(q[l]+[ 1, 1],q[m],v[l],v[m])
            t8=compute_dt(q[l]+[ 1,-1],q[m],v[l],v[m])
            t9=compute_dt(q[l]+[-1,-1],q[m],v[l],v[m])
                
            a=[t1,t2,t3,t4,t5,t6,t7,t8,t9]
            Collisions[m,l]=minimum([t1,t2,t3,t4,t5,t6,t7,t8,t9])
    
            if isinf(Collisions[m,l])
                Fantome[m,l]=0
            else 
                Fantome[m,l]=argmin(a)
            end
            Fantome[m,l]
     
        end
        
        for l=m+1:N
             
            t1=compute_dt(q[l]+[ 0, 0],q[m],v[l],v[m])
            t2=compute_dt(q[l]+[ 1, 0],q[m],v[l],v[m])
            t3=compute_dt(q[l]+[-1, 0],q[m],v[l],v[m])
            t4=compute_dt(q[l]+[ 0, 1],q[m],v[l],v[m])
            t5=compute_dt(q[l]+[ 0,-1],q[m],v[l],v[m])
            t6=compute_dt(q[l]+[-1, 1],q[m],v[l],v[m]) 
            t7=compute_dt(q[l]+[ 1, 1],q[m],v[l],v[m])
            t8=compute_dt(q[l]+[ 1,-1],q[m],v[l],v[m])
            t9=compute_dt(q[l]+[-1,-1],q[m],v[l],v[m])
                
            a=[t1,t2,t3,t4,t5,t6,t7,t8,t9]
            Collisions[m,l] = minimum([t1,t2,t3,t4,t5,t6,t7,t8,t9])
    
            if isinf(Collisions[m,l])
                Fantome[m,l] = 0
            else
                Fantome[m,l] = argmin(a)
            end
    
         end
         
         for l=1:n-1
    
             t1=compute_dt(q[l]+[ 0, 0],q[n],v[l],v[n])
             t2=compute_dt(q[l]+[ 1, 0],q[n],v[l],v[n])
             t3=compute_dt(q[l]+[-1, 0],q[n],v[l],v[n])
             t4=compute_dt(q[l]+[ 0, 1],q[n],v[l],v[n])
             t5=compute_dt(q[l]+[ 0,-1],q[n],v[l],v[n])
             t6=compute_dt(q[l]+[-1, 1],q[n],v[l],v[n]) 
             t7=compute_dt(q[l]+[ 1, 1],q[n],v[l],v[n])
             t8=compute_dt(q[l]+[ 1,-1],q[n],v[l],v[n])
             t9=compute_dt(q[l]+[-1,-1],q[n],v[l],v[n])
                
             a = [t1,t2,t3,t4,t5,t6,t7,t8,t9]
             Collisions[n,l] = minimum([t1,t2,t3,t4,t5,t6,t7,t8,t9])
    
             if isinf(Collisions[n,l])
                 Fantome[n,l] = 0
             else
                 Fantome[n,l] = argmin(a)
             end
    
          end
         
          for l=n+1:N
    
              t1=compute_dt(q[l]+[ 0, 0],q[n],v[l],v[n])
              t2=compute_dt(q[l]+[ 1, 0],q[n],v[l],v[n])
              t3=compute_dt(q[l]+[-1, 0],q[n],v[l],v[n])
              t4=compute_dt(q[l]+[ 0, 1],q[n],v[l],v[n])
              t5=compute_dt(q[l]+[ 0,-1],q[n],v[l],v[n])
              t6=compute_dt(q[l]+[-1, 1],q[n],v[l],v[n]) 
              t7=compute_dt(q[l]+[ 1, 1],q[n],v[l],v[n])
              t8=compute_dt(q[l]+[ 1,-1],q[n],v[l],v[n])
              t9=compute_dt(q[l]+[-1,-1],q[n],v[l],v[n])
                
              a=[t1,t2,t3,t4,t5,t6,t7,t8,t9]
              Collisions[n,l] = minimum([t1,t2,t3,t4,t5,t6,t7,t8,t9])
    
              if isinf(Collisions[n,l])
                  Fantome[n,l] = 0
              else
                  Fantome[n,l] = argmin(a)
              end
    
          end
              
    end

    gif(anim, joinpath(@__DIR__, "hs_periodic.gif"), fps = 10)
end

@time main(100)
