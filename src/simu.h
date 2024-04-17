#pragma once

#include "mesh.h"
#include "fem.h"
#include <math.h>
#include <cmath>
#include <iostream>

namespace FEM2A {
    namespace Simu {

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

        //#################################
        //  Simulations
        //#################################

        void pure_dirichlet_pb( const std::string& mesh_filename, bool verbose )
        {
        	std::cout << "Solving a pure Dirichlet problem" << std::endl;
        	
	    	Mesh M;
		M.load("data/square.mesh");
		
        	ShapeFunctions ShFct_triangle(2, 1);
        	Quadrature quad_triangle = Quadrature::get_quadrature(2,false);
        	
        	SparseMatrix K(M.nb_vertices());
        	
            	std::vector< bool > attribute_is_dirichlet(2,false);
            	attribute_is_dirichlet[1] = true;
            	
            	std::vector< double > values(M.nb_vertices(), 0);
        	M.set_attribute(unit_fct, 1, true);
            
            
        	for (int triangle = 0; triangle < M.nb_triangles(); triangle++) {
        		ElementMapping ElMap(M, false, triangle);
        		DenseMatrix Ke;
        		assemble_elementary_matrix(ElMap, ShFct_triangle, quad_triangle, unit_fct, Ke);
     			local_to_global_matrix(M, triangle, Ke, K);
        	}


        	for (int vertice = 0; vertice < M.nb_vertices(); vertice++) {
 			values[vertice] = xy_fct(M.get_vertex(vertice));
        	}
    
   
           	if ( verbose ) {
                	std::cout << " K matrix \n" << std::endl;
                	K.print();
            	}
            	
            	std::vector< double > F(M.nb_vertices(), 0.);

            	apply_dirichlet_boundary_conditions(M, attribute_is_dirichlet, values,K, F); 
            	
		std::vector<double> x(M.nb_vertices(), 0);
		
		solve(K, F, x);
		
		M.save("Pur_Dirichlet.mesh");
		save_solution( x, "Pur_Dirichlet.bb" ) ;
		
		            	
        }

    }

}
