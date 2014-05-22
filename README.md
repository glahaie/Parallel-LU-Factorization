INF7235_TP2
===========

Approches de parallélisation par échange de message pour la décomposition LU. Il y a
quatre versions disponibles:

 - Version séquentiel
 - Version avec agglomération par lignes adjacentes
 - Version avec agglomération par lignes cycliques
 - Version avec agglomération 2D adjacente


**Utilisation**

Arguments pour tous les programmes:

 - Obligatoire: -m <taille_matrice>: taille de la matrice à traiter.On doit le fournir même
   si on traite un fichier, et les grandeurs doivent être identiques. Le programme ne vérifie pas ça
 - -f <cheminfichier>: fournir un fichier contenant une matrice précise. S'il n'est pas
   fourni, le programme génère une matrice avec des données aléatoires
 - -s <cheminfichier>: fichier de solution, on veut vérifier l'exactitude du calcul
 - -i: impression à l'écran du résultat
 - -w <nomFichier>: écriture du résultat dans le fichier
 - -p: demander au programme de faire l'algorithme avec un pivot partiel. Pas disponible pour lu\_mpi\_cart
 - -x, y <nombre>: pour lu\_mpi\_cart seulement, on doit obligatoirement spécifier le nombre de lignes / colonnes
   de blocs pour le programme. Ces nombres doivent diviser le nombre de noeuds.

Il est possible de tester l'exactitude des différents programmes avec make tests, pour la performance,
on peut faire make mesures NP=<nombre> MATRIX\_SIZE=<nombre>. Pour lu\_mpi\_cart, on doit aussi fournir
ROWS et COLS
