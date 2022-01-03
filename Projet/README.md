
# lscpu
Architecture:                    x86_64
CPU op-mode(s):                  32-bit, 64-bit
Byte Order:                      Little Endian
Address sizes:                   39 bits physical, 48 bits virtual
CPU(s):                          4
On-line CPU(s) list:             0-3
Thread(s) per core:              2
Core(s) per socket:              2
Socket(s):                       1
NUMA node(s):                    1
Vendor ID:                       GenuineIntel
CPU family:                      6
Model:                           61
Model name:                      Intel(R) Core(TM) i5-5300U CPU @ 2.30GHz
Stepping:                        4
CPU MHz:                         881.367
CPU max MHz:                     2900,0000
CPU min MHz:                     500,0000
BogoMIPS:                        4589.52
Virtualization:                  VT-x
L1d cache:                       64 KiB
L1i cache:                       64 KiB
L2 cache:                        512 KiB
L3 cache:                        3 MiB
NUMA node0 CPU(s):               0-3

# Mesure du temps d'exécution en séquentiel

![](https://markdown.data-ensta.fr/uploads/upload_e4e723172824868191ec43a245de2a1d.png)


# Parallélisation affichage contre simulation

## Méthode suivie

Afin d'afficher les résultats de la simulation sans que cette dernière soit retardée, on choisit de paralléliser la partie affichage et la partie simulation à travers MPI.
On a choisi la solution simple qui consiste à faire communiquer les deux processus en envoyant deux std::vector<int> qui contiennent les valeurs nécessaires à afficher. L'extraction de ces valeurs des objets Grille qui se faisait avant dans la fonction d'affichage se fera désormais dans la partie simulation.

## Version synchrone contre la version asynchrone
    
![](https://markdown.data-ensta.fr/uploads/upload_88ad679f525cb1d7e2bdd3a5e3b3507e.png)

----
D'après la méthode choisie, il était attendu que le résultat de la parallélisation ne soit pas aussi rapide que la version séquentiel sans affichage. 
On remarque aussi que la version asynchrone est plus rapide. Ceci revient à une communication moins fréquente entre le processus de simulation et celui de l'affichage. De plus, on a remarqué qu'on affiche les résultats chaque 2 jours en asynchrone (pour 100 000 individus), ce qui est cohérent car le temps d'exécution de la simulation avec l'affichage est inférieur à deux fois le temps d'exécution sans affichage (fig. 1). 

On calcule les speed-ups par rapport à la version en séquentiel (avec affichage) : 
| Population | Speed-up en sync | Speed-up en async |
| ---------- | ---------------- | ----------------- |
| 100 000     | 1.46             | 1.76              |
| 200 000     | 1.48             | 1.63              |
| 400 000     | 1.51             | 1.63              |


# Parallélisation de la simulation avec OpenMP
    
![](https://markdown.data-ensta.fr/uploads/upload_b7fbcfb9d392542840e0fed96dd84049.png)
    
Afin de paralléliser la partie simulation du programme, on utilise Open MP pour paralléliser les boucles dans la fonction *simulation* et *màjStatistique* en protégeant les variables nécessaires avec un *reduce (+: )* . 
On remarque que l'effet de parallélisation avec Open MP est de plus en plus notable quand on choisit des populations plus larges.

On calcule les speed-ups de la version OMP avec 4 threads par rapport à la version en séquentiel (avec affichage):


| Population | Speed-up |
| ---------- | -------- |
| 100 000    | 1.64     |
| 200 000    | 1.85     |
| 400 000    | 1.98     |

Paradoxalement, paralléliser avec Open MP en utilisant 2, 3 ou 4 threads donne le même résultat à peu près. Il est possible que le problème est devenu rapidement memory bound puisque les boucles ne contiennent pas des instructions très complexes. Mais, il se peut aussi qu'il soit un problème matériel puisque l'ordinateur où les essais sont effectués n'est pas très puissant.

![](https://markdown.data-ensta.fr/uploads/upload_210c2a2de06657727bc1219c79f00987.png)
    
On peut calculer également le speed-up par rapport à la simulation initiale en prenant un nombre constant comme population/thread:
    

| Pop / thread | 2 threads | 4 threads |
| ------------ | --------- | --------- |
| 50 000       | 1.59      | 1.85      |
| 100 000      | 1.81      | 1.98      |

On remarque que le speed-up avec une population/thread constante s'améliore si on prend un nombre de threads plus grand et en prenant un nombre d'individus plus grand.
    
# Parallélisation de la partie simulation avec MPI

Après Open MP, on veut paralléliser la partie simulation du programme avec MPI. L'idée consiste à partager les éléments du *std::vector<Individu> population* entre les différents processus s'occupant de la simulation. Chaque processus aura donc un nombre d'individus pour lesquels il doit calculer s'ils sont contaminés ou pas en utilisant la grille qui sera commune (dupliquée puisqu'on est en mémoire distribuée), et puis mettre à jour cette grille.
!! La version disponible de cette parallélisation ne marche pas malheureusement pour plus qu'un seul processeur pour la partie simulation. !!