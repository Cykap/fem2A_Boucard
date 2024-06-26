#pragma once

#include "mesh.h"
#include "fem.h"
#include "solver.h"

#include <assert.h>
#include <iostream>
#include <iomanip>
#include <cmath>
#include <algorithm>
#include <stdlib.h>

namespace FEM2A {
    namespace Tests {
	//#################################
        //  Useful functions
        //#################################

        double unit_fct( vertex v )
        {
            return 1.;
        }

        double zero_fct( vertex v )
        {
            return 0.;
        }

        double xy_fct( vertex v )
        {
            return v.x + v.y;
        }
	
	
        bool test_load_mesh()
        {
            Mesh mesh;
            mesh.load("data/square.mesh");

            std::cout << "Vertices <x> <y> <att>" << std::endl;
            for( int v = 0; v < mesh.nb_vertices(); v++ ) {
                std::cout << mesh.get_vertex(v).x << " " << mesh.get_vertex(v).y
                    << " " << mesh.get_vertex_attribute(v) << std::endl;
            }

            std::cout << "Edges <v0> <v1> <att>" << std::endl ;
            for( int ed = 0; ed < mesh.nb_edges(); ed++ ) {
                std::cout << mesh.get_edge_vertex_index(ed, 0) << " "
                    << mesh.get_edge_vertex_index(ed, 1) << " "
                    << mesh.get_edge_attribute(ed) << std::endl;
            }

            std::cout << "Triangles <v0> <v1> <v2> <att>" << std::endl ;
            for( int tr = 0; tr < mesh.nb_triangles(); tr++ ) {
                std::cout << mesh.get_triangle_vertex_index(tr, 0) << " "
                    << mesh.get_triangle_vertex_index(tr, 1) << " "
                    << mesh.get_triangle_vertex_index(tr, 2) << " "
                    << mesh.get_triangle_attribute(tr) << std::endl;
            }

            return true;
        }

        bool test_load_save_mesh()
        {
            Mesh mesh;
            mesh.load("data/geothermie_4.mesh");
            mesh.save("data/geothermie_4.mesh");
            return true;
        }
        
        bool test_quadrature()
        {
        	Quadrature quad;
        	quad = quad.get_quadrature(2,false);
        	int i=0;
        	float somme = 0;
        	while (i<quad.nb_points()) {
        		somme = somme + quad.weight(i);
        		i = i+1;
        	}
        	std::cout << "somme : " << somme << std::endl;
        	return true;
        }

	bool test_ElementMapping()
	{
		Mesh carre;
		carre.load("data/square.mesh");
		// test sur le triangle 4 de square.mesh
		ElementMapping( carre, false, 4 );
		return true;
	}
	
	bool test_ElementTransform()
	{
		Mesh carre;
		carre.load("data/square.mesh");
		vertex p;
		p.x = 0.2;
		p.y = 0.4;
		vertex r;
		ElementMapping E( carre, false, 4 );
		r = E.transform(p);
		std::cout << r.x << "	" << r.y << std::endl;
		return true;
	}
	
	bool test_JacobMatrix()
	{
		Mesh carre;
		carre.load("data/square.mesh");
		vertex p;
		p.x = 0.2;
		p.y = 0.4;
		DenseMatrix J;
		ElementMapping E( carre, false, 4 );
		J = E.jacobian_matrix(p);
		J.print();
		double det = E.jacobian(p);
		std::cout << det << std::endl;
		return true;
	}
	
	bool test_ElementaryMatrix() 
	{
        	Mesh carre;
		carre.load("data/square.mesh");
		ElementMapping ElMap( carre, false, 4 );
		ShapeFunctions ShFct(2, 1);
		Quadrature quad;
        	quad = quad.get_quadrature(2,false);
        	DenseMatrix Ke;
        	assemble_elementary_matrix(ElMap, ShFct, quad, unit_fct, Ke);
        	Ke.print();
        	SparseMatrix K(carre.nb_vertices());
        	local_to_global_matrix(carre, 4, Ke, K);
        	K.print();
        	
        	return true;
        }
        
        bool test_ElementaryVector() 
	{
        	Mesh carre;
		carre.load("data/square.mesh");
		ElementMapping ElMap( carre, false, 4 );
		ShapeFunctions ShFct(2, 1);
		Quadrature quad;
        	quad = quad.get_quadrature(2,false);
        	std::vector< double > Fe(ShFct.nb_functions(),0);
        	assemble_elementary_vector(ElMap, ShFct, quad, unit_fct, Fe);
		std::cout << Fe[0] << std::endl;
		std::cout << Fe[1] << std::endl;
        	std::cout << Fe[2] << std::endl;
        	return true;
        }
        
        bool test_NeumannVector() 
	{
        	Mesh carre;
		carre.load("data/square.mesh");
		ElementMapping ElMap( carre, true, 4 );
		ShapeFunctions ShFct(1, 1);
		Quadrature quad;
        	quad = quad.get_quadrature(2, true);
        	std::vector< double > Fe(ShFct.nb_functions(),0);
        	assemble_elementary_neumann_vector(ElMap, ShFct, quad, unit_fct, Fe);
		std::cout << Fe[0] << std::endl;
		std::cout << Fe[1] << std::endl;
        	std::cout << Fe[2] << std::endl;
        	return true;
        }
        
    }
}
