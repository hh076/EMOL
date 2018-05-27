#!/usr/bin/ruby

printf "ruby version: %s\n", RUBY_VERSION
if ( RUBY_VERSION.to_f <= 1.8 ) then
  def require_relative( str )
    str2 = "./" + str
    require str2
  end
end

require_relative 'emol_util'
require 'gsl'

class Liu
   def initialize( _a, _nsol, _thresh )
      @a = _a.clone
      @nsol = _nsol
      @thresh = _thresh

#                                           (1) Trial vectors
      @n = @a.get_n
      @b = @a.init_vec( _nsol )
      @ab = @a.get_hc( @b )
#     show_matrix( @b, "initial b-matrix" )
#     show_matrix( @ab, "initial ab-matrix" )

      if EmolConsts::PRINT_LEVEL > 1
      	printf( "\n\n*** Davidson-Liu method Diagonalization ***\n" )  
      	printf( "Number of initial trial vectors : %5d\n", @b.size2 )  
      	printf( "Convergence threshold :  %14.10f \n", @thresh )  
      	printf( "Matrix size :  %10d \n", @n )  
      	printf( "Number of roots :  %10d \n", @nsol )  
      end
   end

   def solve
      if EmolConsts::PRINT_LEVEL > 1
      	printf( " \n" )  
      	printf( "     n    dim    E                 Elow             Dev \n" )  
      end	
      iteration = -1; m = 0
      prev_e = Array.new( @nsol, 0.0 )

      n_add = @b.size2
      while n_add > 0
	 iteration += 1
#                                           (2) Small matrix eigenvalue problem
         m = @b.size2
         rdmat = GSL::Matrix.alloc(m, m)
         for i in 0...m
            for j in 0..i
	       rdmat[i, j] = @b.col(i).col * @ab.col(j)
	    end
         end 
#        EmolUtil::print_2dim_ary( rdmat.to_a, "Small Matrix", 8, "%12.6f" )
         eval, evec = GSL::Eigen::symmv( rdmat )
         GSL::Eigen::symmv_sort( eval, evec, type=GSL::Eigen::SORT_VAL_ASC )
#        EmolUtil::print_2dim_ary_with_value( evec.to_a, eval.to_a, "Eigenvalue and Vectors of Small Matrix", 8, "%12.6f" )

#                                           (3) Correction vectors
         f = []
         for k in 0...@nsol
	    d = GSL::Vector.alloc( @n )
	    vec = evec.col( k ); val = eval[ k ]
            for i in 0...m
               d += vec[ i ] * ( @ab.col( i ) - val * @b.col( i ) )
	    end

	    for p in 0...@n
               d[ p ] = d[ p ] / ( val - @a.get_diag( p ) )    
	    end
	    f.push( d )
         end
#         show_vectors( f, "Correction vectors" )

#                                           Print energy lowerings
         for k in 0...@nsol
	    lowering = eval[ k ] - prev_e[ k ]
#     p f[ k ] * f[ k ].col 
      	    if EmolConsts::PRINT_LEVEL > 1
           	printf( " %5d  %5d  %14.10f    %14.10f     %10.8f \n", 
	   	iteration, m, eval[ k ], lowering, Math::sqrt( f[ k ] * f[ k ].col ) )
	    end
	    prev_e[ k ] = eval[ k ]
	 end

#                                           (4, 5) Check convergence and Schmidt orthonormaization 
	 n_add = 0
         for k in 0...@nsol
            @b.each_col do | b |
               f[ k ] = f[ k ] - ( f[ k ] * b ) * b.col
	    end
	    if Math::sqrt( f[ k ] * f[ k ].col ) > @thresh
	       @b = @b.horzcat( GSL::Matrix.alloc( f[ k ].normalize.col ) )
	       n_add += 1
	    end
         end
#        show_matrix( @b, "New element vectors" )

#                                           (6) Ab generation for newly added b vectors
         if n_add > 0
	    @ab = @ab.horzcat( @a.get_hc( @b.submatrix( nil, @b.size2 - n_add...@b.size2 ) ) )
	 end
      end

      return eval, @b * evec.submatrix( nil, 0...@nsol )
   end

   def show_vectors( v, title )
      print "\n", title, "\n"
      for p in 0...@n do
         printf( "%5d ", p + 1 )
	 for k in 0...v.size do
            printf( "%12.6f ", v[k][p] )
	 end
         print "\n"
      end
   end

   def show_matrix( v, title )
      print "\n", title, "\n"
      for p in 0...v.size1 do
         printf( "%5d ", p + 1 )
	 for k in 0...v.size2 do
            printf( "%12.6f ", v[p, k] )
	 end
         print "\n"
      end
   end
end
