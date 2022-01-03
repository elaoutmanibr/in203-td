Population = [1000,2000,5000,10000,20000,50000,100000,200000,400000]
T1 = [3,4,5,7,13,28,51,98,200] # code original avec affichage
T2 =[1,1,1,2,4,11,24,52,114] # code orig sans affichage
T3 =[0,0,0,0,0,17,35,66,132] # parallelisation MPI sync
T4 = [0,0,0,0,0,14,29,60,123] # parallelisation MPI async
T5 =[0,0,0,0,0,18,31,53,101] # OMP 4


import matplotlib.pyplot as plt

plt.plot(Population,T1,'r')
plt.plot(Population,T2,'b')
plt.plot(Population,T3,'g')
plt.plot(Population,T4,'k')
plt.plot(Population,T5,'m')
plt.title("Temps d'exécution des différents programmes")
plt.xlabel("Population")
plt.ylabel("Temps d'exécution en ms")
plt.legend(["Séquentiel","Séquentiel sans affichage","Parallélisation MPI (sync affichage)","Parallélisation MPI (async affichage)","Parallélisation MPI async + OMP (4 threads)"], loc = "lower right")
plt.show()

Population = [50000,100000,200000,400000,800000]
S1 = [18,32,54,104,202] # OMP 2
S2 = [18,29,54,100,202] # OMP 3
S3 = [17,29,53,101,202] # OMP 4
plt.plot(Population,S1,'r')
plt.plot(Population,S2,'b')
plt.plot(Population,S3,'g')

plt.legend(["2 threads","3 threads","4 threads"],loc="lower right")
plt.show()
