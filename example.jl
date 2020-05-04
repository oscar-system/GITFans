############################################################################
##  The computations shown here belong to Example 5.2 in the paper at
##  https://arxiv.org/abs/1603.09241
##  (Computing GIT-Fans with Symmetry and the Mori Chamber Decomposition ...)

using GITFans
using GAP
using Polymake
using Singular


############################################################################
##  1. the input data (see the paper)

# grading matrix
Q = [
 1  1   0   0   0 ;
 1  0   1   1   0 ;
 1  0   1   0   1 ;
 1  0   0   1   1 ;
 0  1   0   0  -1 ;
 0  1   0  -1   0 ;
 0  1  -1   0   0 ;
 0  0   1   0   0 ;
 0  0   0   1   0 ;
 0  0   0   0   1 ];

# polynomial ring
R, T = Singular.PolynomialRing( Singular.QQ,
           [ "x" * string(i) for i in 1:size( Q, 1 ) ] )

# generators for the ideal
II = Singular.Ideal( R, [
    T[5]*T[10] - T[6]*T[9] + T[7]*T[8],
    T[1]*T[9]  - T[2]*T[7] + T[4]*T[5],
    T[1]*T[8]  - T[2]*T[6] + T[3]*T[5],
    T[1]*T[10] - T[3]*T[7] + T[4]*T[6],
    T[2]*T[10] - T[3]*T[9] + T[4]*T[8],
] )

# symmetries: (2,3)(5,6)(9,10), (1,5,9,10,3)(2,7,8,4,6)
perms_list = [ [1,3,2,4,6,5,7,8,10,9], [5,7,1,6,9,2,8,4,10,3] ];
perms = [ GAP.Globals.PermList(GAP.julia_to_gap(i)) for i in perms_list ];
G = GAP.Globals.Group( GAP.julia_to_gap( perms ) )


############################################################################
##  2. preprocessing: compute G-orbits and induced action

fan_descr = GITFans.GITFan( II, Q, G );


############################################################################
##  3. fan traversal: compute neighbors until the list is complete

hash_list = GITFans.fan_traversal( fan_descr, Q );

length( hash_list )


############################################################################
##  4. translate the hashes back to cones

orbit_list = fan_descr[1];

result_cones = map(x -> Polymake.polytope.intersection(
                          GITFans.cones_from_bitlist(orbit_list, x)...), hash_list);


############################################################################
##  5. expand the orbits (yields 76 maximal cones)

hom = GITFans.action_on_target( Q, G );
expanded = GITFans.orbit_cone_orbits( result_cones, hom; disjoint_orbits = false )
length( expanded )              # 6
map( length, expanded )
sum( map( length, expanded ) )  # 76


############################################################################
##  6. create the intersection graph

maxcones = [];
for orb in expanded
  append!( maxcones, orb )
end

dims = Set( map( Polymake.polytope.dim, maxcones ) )
codim1 = collect( dims )[1] - 1

edges = []
for i in 1:length( maxcones )
  for j in 1:(i-1)
    if Polymake.polytope.dim(
           Polymake.polytope.intersection( maxcones[i], maxcones[j] ) ) == codim1
      push!( edges, [ j, i ] )
    end
  end
end

length( edges )    # 180

edges = convert( Vector{Vector{Int}}, edges );

intergraph = Polymake.graph.graph_from_edges( edges );

Polymake.graph.visual( intergraph )

