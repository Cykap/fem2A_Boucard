Ce document a pour objectif de décrire comment réaliser les tests et simulations du projet. Pour se faire, il suffira de passer le bon booléen sur True dans le fichier main.cpp pour effectuer le test ou la simulation souhaité.

Pour les tests :
	t_opennl, t_lmesh et t_io sont les variables déjà présentes sur ce document, elles ne sont pas utilisées.
	t_quadra : Permet de faire le test de quadrature
	t_map : Permet de faire le test de ElementMapping (déterminer les coordonnées des sommets du triangle)
		-> Nécessite de réactiver l'affichage dans la boite de dialogue dans la fonction ElementMapping du fichier fem.cpp
	t_transfo : Permet de faire le test de transform de ElementMapping
	t_jacobM : Permet de faire le test de la matrice jacobienne d'un point et son déterminant
	t_elemMatrix : Permet de tester la création d'une matrice élémentaire
	t_elemVect : Permet de tester la création d'un vecteur élémentaire 
	t_neumann : Permet de tester la création d'un vecteur élémentaire de Neumann

Pour les simulations :
	simu_pure_dirichlet : Permet de simuler le problème de Dirichlet pur pour le maillage carre.mesh
	simu_dirichlet_source : Permet de simuler le problème de Dirichlet avec un terme source pour le maillage carre.mesh
	simu_dirichlet_sin : Permet de simuler le problème sinus bump pour le maillage carre.mesh
	simu_dirichlet_neum : Permet de simuler le problème de Dirichlet et Neumann pour le maillage carre.mesh
	simu_mug : Permet de simuler le problème du mug pour le maillage mug_1.mesh
	simu_geo : Permet de simuler le problème de géothermie pour le maillage geothermie_0_1.mesh

	
