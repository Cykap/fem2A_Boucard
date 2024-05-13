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
        
        double hundred_fct( vertex v )
        {
            return 100.;
        }
        
        double decimal_fct( vertex v )
        {
            return 0.1;
        }

        double xy_fct( vertex v )
        {
            return v.x + v.y;
        }

        double sinus_fct( vertex v )
        {
            double PI = 3.1415;
            return 2 * PI * PI * sin(PI*v.x) * sin(PI*v.y);
        }
        
        double sinus_sol_fct( vertex v )
        {
            double PI = 3.1415;
            return sin(PI*v.x) * sin(PI*v.y);
        }
        
        double neumann_fct( vertex v )
        {
            double PI = 3.14159265;
            return sin(PI*v.y);
        }
        
        double test_droite( vertex v )
        {
            if (abs(v.x - 1)<0.00001) {
            	return 1.;
            }
            else {
            	return 0.;
            }
        }
        
        double test_gauche( vertex v )
        {
            if (abs(v.x - 0)<0.00001) {
            	return 1.;
            }
            else {
            	return 0.;
            }
        }
        
        double test_mug_dirichlet( vertex v )
        {
            if (abs(v.x - 1)<0.00001) {
            	if ( (v.y >= 1) && (v.y <= 10)) {
            		return 1.;
            	}

            }
            else if (abs(v.x - 20)<0.00001) {
		if ( (v.y >= 1) && (v.y <= 10)) {
            		return 1.;
            	}
            }
            else if (abs(v.y - 1)<0.00001) {
		if ( (v.x >= 1) && (v.x <= 20)) {
            		return 1.;
            	}
            }
            else {
            	return 0.;
            }
            
       	return 0.;
        }
        
        double neumann_mug(vertex v)
        {
        	return -0.1;
        }
        
	double test_geo_dirichlet( vertex v )
        {
            if (abs(v.y - 20)<0.00001) {
            		return 1.;
            }

            else {
            	return 0.;
            }
            
        }
        

        //#################################
        //  Simulations
        //#################################

        void pure_dirichlet_pb( const std::string& mesh_filename, bool verbose )
        {
        	std::cout << "Solving a pure Dirichlet problem" << std::endl;
        	
	    	Mesh M;
		M.load(mesh_filename);
		
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
		
		M.save("Pur_Dirichlet_fine_mesh.mesh");
		save_solution( x, "Pur_Dirichlet_fine_mesh.bb" ) ;           	
        }
        
        
        void dirichlet_source_pb( const std::string& mesh_filename, bool verbose )
        {
        	std::cout << "Solving a Dirichlet problem with source" << std::endl;
        	
	    	Mesh M;
		M.load(mesh_filename);
		
        	ShapeFunctions ShFct_triangle(2, 1);
        	Quadrature quad_triangle = Quadrature::get_quadrature(2,false);
        	
        	SparseMatrix K(M.nb_vertices());
        	std::vector< double > F(M.nb_vertices(), 0.);
        	
            	std::vector< bool > attribute_is_dirichlet(2,false);
            	attribute_is_dirichlet[1] = true;
            	
            	std::vector< double > values(M.nb_vertices(), 0);
        	M.set_attribute(unit_fct, 1, true);
            
            
        	for (int triangle = 0; triangle < M.nb_triangles(); triangle++) {
        		ElementMapping ElMap(M, false, triangle);
        		DenseMatrix Ke;
        		assemble_elementary_matrix(ElMap, ShFct_triangle, quad_triangle, unit_fct, Ke);
     			local_to_global_matrix(M, triangle, Ke, K);
     			
     			std::vector< double > Fe(ShFct_triangle.nb_functions(),0);
     		        assemble_elementary_vector(ElMap, ShFct_triangle, quad_triangle, unit_fct, Fe);
     		        local_to_global_vector(M, false, triangle, Fe, F);
        	}


        	for (int vertice = 0; vertice < M.nb_vertices(); vertice++) {
 			values[vertice] = zero_fct(M.get_vertex(vertice));
        	}
    
   
           	if ( verbose ) {
                	std::cout << " K matrix \n" << std::endl;
                	K.print();
            	}

            	apply_dirichlet_boundary_conditions(M, attribute_is_dirichlet, values, K, F); 
            	
		std::vector<double> x(M.nb_vertices(), 0);
		
		solve(K, F, x);
		
		M.save("Dirichlet_source_fine_mesh.mesh");
		save_solution( x, "Dirichlet_source_fine_mesh.bb" ) ;           	
        }
        
        
        void dirichlet_sinus_pb( const std::string& mesh_filename, bool verbose )
        {
        	std::cout << "Solving a Dirichlet problem with sinus bump" << std::endl;
        	
	    	Mesh M;
		M.load(mesh_filename);
		
        	ShapeFunctions ShFct_triangle(2, 1);
        	Quadrature quad_triangle = Quadrature::get_quadrature(2,false);
        	
        	SparseMatrix K(M.nb_vertices());
        	std::vector< double > F(M.nb_vertices(), 0.);
        	
            	std::vector< bool > attribute_is_dirichlet(2,false);
            	attribute_is_dirichlet[1] = true;
            	
            	std::vector< double > values(M.nb_vertices(), 0.);
        	M.set_attribute(unit_fct, 1, true);
            
            
        	for (int triangle = 0; triangle < M.nb_triangles(); triangle++) {
        		ElementMapping ElMap(M, false, triangle);
        		DenseMatrix Ke;
        		assemble_elementary_matrix(ElMap, ShFct_triangle, quad_triangle, unit_fct, Ke);
     			local_to_global_matrix(M, triangle, Ke, K);
     			
     			std::vector< double > Fe(ShFct_triangle.nb_functions(),0.);
     		        assemble_elementary_vector(ElMap, ShFct_triangle, quad_triangle, sinus_fct, Fe);
     		        local_to_global_vector(M, false, triangle, Fe, F);
        	}


        	for (int vertice = 0; vertice < M.nb_vertices(); vertice++) {
 			values[vertice] = zero_fct(M.get_vertex(vertice));
        	}
    
   
           	if ( verbose ) {
                	std::cout << " K matrix \n" << std::endl;
                	K.print();
            	}

            	apply_dirichlet_boundary_conditions(M, attribute_is_dirichlet, values, K, F); 
            	
		std::vector<double> x(M.nb_vertices(), 0);
		
		std::vector<double> x_sol(M.nb_vertices(), 0);
		
		solve(K, F, x);
		
		bool sol = false;
		if (sol == true) {
			for (int i = 0; i < M.nb_vertices(); ++i){
				vertex v = M.get_vertex(i);
				x[i] -= sinus_sol_fct( v );
			}
		}
		
		M.save("Dirichlet_sinus.mesh");
		save_solution( x, "Dirichlet_sinus.bb" ) ;           	
        }
        
        
        void dirichlet_neumann_pb( const std::string& mesh_filename, bool verbose )
        {
        	std::cout << "Solving a Dirichlet problem with Neumann" << std::endl;
        	
	    	Mesh M;
		M.load(mesh_filename);
		
        	ShapeFunctions ShFct_triangle(2, 1);
        	ShapeFunctions ShFct_1D(1, 1);
        	Quadrature quad_triangle = Quadrature::get_quadrature(2,false);
        	Quadrature quad_1D = Quadrature::get_quadrature(2,true);
        	
        	SparseMatrix K(M.nb_vertices());
        	std::vector< double > F(M.nb_vertices(), 0.);
        	
            	std::vector< bool > attribute_is_dirichlet(2,false); 
            	attribute_is_dirichlet[1] = true;
            	
            	std::vector< double > values(M.nb_vertices(), 0.); // 0 : None, 1 : Dirichlet, 2 : Neumann
            	M.set_attribute(unit_fct, 0, true);
        	M.set_attribute(test_droite, 1, true);
        	M.set_attribute(test_gauche, 2, true);
          
        	for (int triangle = 0; triangle < M.nb_triangles(); triangle++) {
        		ElementMapping ElMap(M, false, triangle);
        		DenseMatrix Ke;
        		assemble_elementary_matrix(ElMap, ShFct_triangle, quad_triangle, unit_fct, Ke);
     			local_to_global_matrix(M, triangle, Ke, K);
     			
     			std::vector< double > Fe(ShFct_triangle.nb_functions(),0.);
     		        assemble_elementary_vector(ElMap, ShFct_triangle, quad_triangle, unit_fct, Fe);
     		        local_to_global_vector(M, false, triangle, Fe, F);
        	}

        	for (int vertice = 0; vertice < M.nb_vertices(); vertice++) {
 			values[vertice] = zero_fct(M.get_vertex(vertice));
        	}
    
            	apply_dirichlet_boundary_conditions(M, attribute_is_dirichlet, values, K, F); 
    		
		for (int edge; edge < M.nb_edges(); edge++) {
			if (M.get_edge_attribute(edge) == 2) {
			        ElementMapping ElMap(M, true, edge);
				std::vector< double > Fe(ShFct_1D.nb_functions(),0.);
     		        	assemble_elementary_neumann_vector(ElMap, ShFct_1D, quad_1D, neumann_fct, Fe);
     		        	local_to_global_vector(M, true, edge, Fe, F);
			}
		}

		std::vector<double> x(M.nb_vertices(), 0);
		
		solve(K, F, x);
		
		M.save("Dirichlet_neumann.mesh");
		save_solution( x, "Dirichlet_neumann.bb" ) ;
	}   
	
	
	void mug_pb( const std::string& mesh_filename, bool verbose )
        {
        	std::cout << "Solving the mug problem" << std::endl;
        	
	    	Mesh M;
		M.load(mesh_filename);
		
        	ShapeFunctions ShFct_triangle(2, 1);
        	ShapeFunctions ShFct_1D(1, 1);
        	Quadrature quad_triangle = Quadrature::get_quadrature(2,false);
        	Quadrature quad_1D = Quadrature::get_quadrature(2,true);
        	
        	SparseMatrix K(M.nb_vertices());
        	std::vector< double > F(M.nb_vertices(), 0.);
        	
            	std::vector< bool > attribute_is_dirichlet(2,false); 
            	attribute_is_dirichlet[1] = true;
            	
            	std::vector< double > values(M.nb_vertices(), 0.); // 0 : None, 1 : Dirichlet, 2 : Neumann
            	M.set_attribute(unit_fct, 2, true);
        	M.set_attribute(test_mug_dirichlet, 1, true);
          	
        	for (int triangle = 0; triangle < M.nb_triangles(); triangle++) {
        		ElementMapping ElMap(M, false, triangle);
        		DenseMatrix Ke;
        		assemble_elementary_matrix(ElMap, ShFct_triangle, quad_triangle, unit_fct, Ke);
     			local_to_global_matrix(M, triangle, Ke, K);
     		
        	}

        	for (int vertice = 0; vertice < M.nb_vertices(); vertice++) {
        		if (M.get_vertex_attribute(vertice) == 1) {
 				values[vertice] = hundred_fct(M.get_vertex(vertice));
 			}
        	}
    
            	apply_dirichlet_boundary_conditions(M, attribute_is_dirichlet, values, K, F); 
    		
		for (int edge; edge < M.nb_edges(); edge++) {
			if (M.get_edge_attribute(edge) == 2) {
			        ElementMapping ElMap(M, true, edge);
				std::vector< double > Fe(ShFct_1D.nb_functions(),0.);
     		        	assemble_elementary_neumann_vector(ElMap, ShFct_1D, quad_1D, neumann_mug, Fe);
     		        	local_to_global_vector(M, true, edge, Fe, F);
			}
		}

		std::vector<double> x(M.nb_vertices(), 0);
		
		solve(K, F, x);
		
		M.save("Mug_1.mesh");
		save_solution( x, "Mug_1.bb" ) ;
	} 

	void geothermie_pb( const std::string& mesh_filename, bool verbose )
        {
        	std::cout << "Solving the geothermy problem" << std::endl;
        	
	    	Mesh M;
		M.load(mesh_filename);
		
        	ShapeFunctions ShFct_triangle(2, 1);
        	ShapeFunctions ShFct_1D(1, 1);
        	Quadrature quad_triangle = Quadrature::get_quadrature(2,false);
        	Quadrature quad_1D = Quadrature::get_quadrature(2,true);
        	
        	SparseMatrix K(M.nb_vertices());
        	std::vector< double > F(M.nb_vertices(), 0.);
        	
            	std::vector< bool > attribute_is_dirichlet(2,false); 
            	attribute_is_dirichlet[1] = true;
            	
            	std::vector< double > values(M.nb_vertices(), 0.); // 0 : None, 1 : Dirichlet, 2 : Neumann
            	M.set_attribute(unit_fct, 2, true);
        	M.set_attribute(test_geo_dirichlet, 1, true);
          	
          	

        	for (int triangle = 0; triangle < M.nb_triangles(); triangle++) {
        		vertex t1 = M.get_triangle_vertex(triangle, 0);
        		vertex t2 = M.get_triangle_vertex(triangle, 1);
        		vertex t3 = M.get_triangle_vertex(triangle, 2);
        		ElementMapping ElMap(M, false, triangle);
        		DenseMatrix Ke;
        		
        		if ((t1.y <= 5) && (t2.y <= 5) && (t3.y <= 5)) { 
        			assemble_elementary_matrix(ElMap, ShFct_triangle, quad_triangle, unit_fct, Ke);
     				local_to_global_matrix(M, triangle, Ke, K);
     			
     				std::vector< double > Fe(ShFct_triangle.nb_functions(),0.);
     		        	assemble_elementary_vector(ElMap, ShFct_triangle, quad_triangle, hundred_fct, Fe);
     		        	local_to_global_vector(M, false, triangle, Fe, F);
        		}
        		
        		else {
        			assemble_elementary_matrix(ElMap, ShFct_triangle, quad_triangle, decimal_fct, Ke);
     				local_to_global_matrix(M, triangle, Ke, K);
        		}
        	}

        	for (int vertice = 0; vertice < M.nb_vertices(); vertice++) {
        		if (M.get_vertex_attribute(vertice) == 1) {
 				values[vertice] = hundred_fct(M.get_vertex(vertice));
 			}
        	}
    
            	apply_dirichlet_boundary_conditions(M, attribute_is_dirichlet, values, K, F); 
    		
		for (int edge; edge < M.nb_edges(); edge++) {
			if (M.get_edge_attribute(edge) == 2) {
			        ElementMapping ElMap(M, true, edge);
				std::vector< double > Fe(ShFct_1D.nb_functions(),0.);
     		        	assemble_elementary_neumann_vector(ElMap, ShFct_1D, quad_1D, neumann_mug, Fe);
     		        	local_to_global_vector(M, true, edge, Fe, F);
			}
		}

		std::vector<double> x(M.nb_vertices(), 0);
		
		solve(K, F, x);
		
		M.save("Geothermie_0_1.mesh");
		save_solution( x, "Geothermie_0_1.bb" ) ;
	}
        
    }

}
