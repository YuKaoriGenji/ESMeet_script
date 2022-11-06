import numpy as np
import time
from turtle import bgcolor
import xlwt


def create_file_DW(content):
    style_head=xlwt.XFStyle()
    font=xlwt.Font()

    excel=xlwt.Workbook()
    sheet=excel.add_sheet("Quinary Colliding Data")

    head=["k","Exhaustive Searching","Odlyzko's Method","MeetLWE-Rep 0","ES-CS","ES-TS"]
    for index,value in enumerate(head):
        sheet.write(int(0),index,value,style_head)

    for index,value_list in enumerate(content,1):
        for i,value in enumerate(value_list):
            sheet.write(index,i,value)
    file_name=time.time()
    excel.save("./Colliding_hybrid_%s.xls" % file_name)
    return file_name

#security level function
sl(x)=log(x,2)
t(w1,w2)=w1/2+w2

C(a,k)=factorial(a)/(factorial(k)*factorial(a-k))
C2(a,k1,k2)=C(a,k1)*C(a-k1,k2)
C4(a,k1,k2,k3,k4)=C(a,k1)*C(a-k1,k2)*C(a-k1-k2,k3)*C(a-k1-k2-k3,k4)

C_ES(n,q,w1,w2)=C4(n,w1,w1,w2,w2)

test1(n,w1,w2)=C2(n,t(w1,w2),t(w1,w2))
test2(n,w1,w2)=C2(n,t(w1/2,w2/2),t(w1/2,w2/2))
test3(n,w1,w2)=C2(n,t(w1/4,w2/4),t(w1/4,w2/4))
C_NP(d)=(d)^2/(2^1.06)

#================= p_1 =======================
def p1(gsb_norm,m,r,k1,k2):
    s=1
    #print("0 probability: ",(m-2*k1-2*k2)/m)
    for i in range(m-r):
        y1=max_symbolic((gsb_norm[i]-1)/gsb_norm[i],0)
        y2=max_symbolic((gsb_norm[i]-2)/gsb_norm[i],0)
        p_part=(2*(k1/m)*y1+2*(k2/m)*y2+(m-2*k1-2*k2)/m)
        s*=p_part
        #if i%20==0:
            #print("P-Prime Processing ",i," round(s)! \t p_part now: ",p_part.n())
    return s

#================= p_2 =======================
def p2(gsb_norm,m,r,t,k1,k2):
    s=1
    #print("0 probability: ",(m-2*k1-2*k2)/m)
    for i in range(t):
        y1=max_symbolic((gsb_norm[i]-1)/gsb_norm[i],0)
        p_part=2*((k1/2+k2)/m)*y1+(m-2*k1-2*k2)/m
        s*=p_part
        #if i%20==0:
            #print("\tP-Prime Processing ",i," round(s)! \t p_part now: ",p_part.n())
    return s

#================= p_0 =======================
def p0(gsb_norm,m,r,k1,k2):
    s=1
    #print("0 probability: ",(m-2*k1-2*k2)/m)
    #print("m-r: ",m-r)
    for i in range(m-r):
        #y1=max_symbolic((gsb_norm[i]-1)/gsb_norm[i],0)
        #y2=max_symbolic((gsb_norm[i]-2)/gsb_norm[i],0)
        p_part=0
        if 2< gsb_norm[i].n():
            p_part+=2*k2/m
        if 1< gsb_norm[i].n():
            p_part+=2*k1/m

        p_part+=(m-2*k2-2*k1)/m

        s*=p_part
        #if i%20==0:
            #print("\tP-Prime Processing ",i," round(s)! \t p_part now: ",p_part.n())
    return s

def reduction_cost(n,q,k,gamma):
    m=2*n
    d=m
    k2=round(k*gamma/4)*4
    k1=k-k2
    print("k1:",k1,"\tk2:",k2)
    norm_v=(k1+k2*4)^(1/2)
    beta1=0
    print("norm_v:",norm_v.n())
    #------------ calculating the best beta ----------------
    for beta in range(320,500,2):
        cal_v=(d/beta)^(1/2)*((beta*(beta*pi)^(1/beta))/(2*pi*e))^((2*beta-d)/(2*(beta-1)))*(q^n)^(1/d)
        print("[*]beta:",beta,"\tcal_v:",cal_v.n())
        if cal_v>norm_v:
            print("successfully find the beta!")
            beta1=beta
            print("beta1:",beta1)
            break
    Total_T_svp=2^(0.187281*beta1*log(beta1,2)-1.0192*beta1+16.1)
    Total_T_red_under=(d-beta+1)*Total_T_svp
    reduction_cost=Total_T_red_under
    return reduction_cost

def hybrid_cost(n,q,beta,r,target,k,gamma):
    print("================= new ===========================")
    print("r=",r,"\t beta=",beta)
    #n=607
    m=2*n
    #q=18749
    d=m-r
    a=n
    b=n-r
    k2=round(k*gamma/4)*4
    k1=k-k2
    w2=round(k2*r/(m*4))*4
    w1=round(k1*r/(m*4))*4
    print("w1",w1)
    print("w2",w2)
    is_pass=0
    #w1=52-w2
    #target=5
    determinant=1
    #print("About RHF")
    delta=((beta*(beta*pi)^(1/beta))/(2*pi*e))^(1/(2*(beta-1)))
    #print("delta= ",delta.n())
    #k=min_symbolic(floor((b/(log(delta,q)))^(1/2)),d)
    if floor((b/(log(delta,q)))^(1/2))<d:
        k=round(floor((b/(log(delta,q)))^(1/2)).n())
    else:
        k=d
    #print("k= ",k.n())
    gsb_norm=[]
    for ite in range(d):
        if (ite+1)<=d-k:
            gsb_norm.append(q/2)
        else:
            gsb_norm.append(delta^(-2*(ite+1-(d-k)-1)+k)*q^((k-b)/(k))/2)
        determinant*=(2*gsb_norm[ite]).n()
        #if ite%20==0:
            #print("gsb_norm[ ",ite,"]/2= ",gsb_norm[ite].n())
    #print("[**]determinant calculated=",determinant.n())

    p_1=p1(gsb_norm,m,r,k1,k2)
    #print("p_1: ",sl(p_1).n())
    p_2=p2(gsb_norm,m,r,target,k1,k2)
    #print("p_2: ",sl(p_2).n())
    p_0=p0(gsb_norm,m,r,k1,k2)
    #print("p_0: ",sl(p_0).n())
    p_s=(C4(r,w1,w1,w2,w2)*C4(m-r,k1-w1,k1-w1,k2-w2,k2-w2))/(C4(m,k1,k1,k2,k2))*p_0
    #print("for r: ",sl(C4(r,w1,w1,w2,w2)*C4(m-r,k1-w1,k1-w1,k2-w2,k2-w2)).n())
    #print("for m: ",sl(C4(m,k1,k1,k2,k2)).n())
    #print("p_s: ",sl(p_s).n())
    #=============== Reduction Complexity =============
    #print("r:",r,"\tdelta:",delta.n())
    T_svp=2^(0.187281*beta*log(beta,2)-1.0192*beta+16.1)
    T_red_under=(d-beta+1)*T_svp
    #print("T_svp:",log(T_svp,2).n())
    #=============== Colliding Complexity =================
    C_2_list=C4(r,w1/2,w1/4,w2/4,w2/4)
    #print("C_2_list: ",sl(C_2_list).n())
    S_1=C4(r,w1/2,w1/2,w2/2,w2/2)
    R_1=C(w2,w2/2)^(2)*C(w1,w1/2)^(2)
    #print("S_1: ",sl(S_1).n(),"\t R_1: ",sl(R_1).n(),"\t S_1/R_1: ",sl(S_1/R_1).n())
    #print("p_1: ",sl(p_1).n())
    down=1
    for i in range(target):
        down*=2*gsb_norm[i]

    true_S_1=(S_1*3^target)/down
    #print("down: ",sl(down).n())
    #print("true_S_1: ",sl(true_S_1).n())
    #print("true_S_1/R_1: ",sl(true_S_1/R_1).n())
    #---------------- engineering part------------
    C_1=(true_S_1/(p_1*(R_1)))
    ori_C_1=(S_1/(p_1*R_1))
    #---------------- Ori Hybrid Attack-------
    p_p_o=R_1/S_1
    ori_Comp=(S_1/(p_p_o*p_1))^(1/2)+T_red_under
    #print("ori_Comp:",sl(ori_Comp).n())
    #---------------------------------------------------------
    print("[*] criteria:")
    cri_right=(p_2*(C_2_list^2)*(3^target))/(q^target)
    C_1_cri=(true_S_1/(p_1*R_1))
    cri_left=C_1_cri
    print("left=",sl(cri_left).n(),"\t right=",sl(cri_right).n(),"\t We want left smaller")
    is_pass=0
    if(sl(cri_left).n()<=sl(cri_right).n()):
        print("[*]passed!")
        is_pass=1;
    else:
        print("[*]REJECTED!!!!!")
    C_total=2*C_1+2*C_2_list
    #print("C_NP(d):",sl(C_NP(d)).n())
    ES_Comp=(C_total*C_NP(d)/(p_s))+T_red_under
    Hyb_Comp=C_NP(d)*S_1/(p_1*R_1)^(1/2)+T_red_under
    if is_pass==1:
        print("----------------------------")
        print("T_red_under:",log(T_red_under,2).n())
        print("Meet in the Middle:",sl(C_NP(d)*S_1/(p_1*R_1)^(1/2)).n())
        print("C_total:",sl(C_total).n())
        print("----------------------------")
        return (Hyb_Comp,ES_Comp)

    else:
        return (Hyb_Comp,1)


n=607
q=18749
r=210
beta=270
k=200
gamma=0.25
target=8

reduction_result=reduction_cost(n,q,k,gamma)

#for i in range(40,400,40):
    #result=hybrid_cost(n,q,beta,r,target,i,gamma)
    #print("k:",i,"\t ES_Comp:",sl(result).n())

result=hybrid_cost(n,q,beta,r,target,k,gamma)
print("[*]ES_Comp:",sl(result[1]).n())
print("[*]Hyb_Comp:",sl(result[0]).n())
ES_all_result=np.zeros((10,15),dtype=float)
Hyb_all_result=np.zeros((10,15),dtype=float)
for i in range(10):
    for j in range(15):
        now_result=hybrid_cost(n,q,(i+1)*40,(j+1)*20,target,k,gamma)
        ES_all_result[i][j]=sl(now_result[1]).n();
        Hyb_all_result[i][j]=sl(now_result[0]).n();
        print("beta=",(i+1)*40,"r=",(j+1)*20,"\tES_Comp=",ES_all_result[i][j],"\tHyb_Comp=",Hyb_all_result[i][j])
        print("reduction result:",sl(reduction_result).n())

for i in range(10):
    for j in range(15):
        print("beta=",(i+1)*40,"r=",(j+1)*20,"\tES_Comp=",ES_all_result[i][j],"\tHyb_Comp=",Hyb_all_result[i][j])

print("reduction result:",sl(reduction_result).n())

