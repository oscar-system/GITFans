############################################################################
##  The computations shown here belong to Example 5.2 in the following paper:
##  [J. Boehm, S. Keicher, Y. Ren:
##  Computing GIT-Fans with Symmetry and the Mori Chamber Decomposition of
##  $\bar{M}_{0,6}$](https://arxiv.org/abs/1603.09241)

using GITFans
using GAP
using Polymake
using Singular


############################################################################
##  Enter the input as described in Example 5.2 of the paper.

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
G = GAP.Globals.Group(GAP.julia_to_gap(perms))


############################################################################
##  Preprocessing:
##  Encode the GIT fan by the array of orbit cone orbits,
##  the array of homomorphism objects describing the induced actions
##  of G on the orbits,
##  and the matrix Q.

fan_descr = GITFans.GITFan(II, Q, G);


############################################################################
##  Execute the fan traversal,
##  that is, compute neighbors of the known cones,
##  and then the smallest representatives in the G-orbits;
##  each cone is encoded by the set of (indices of) maximal cones
##  containing it.
##  At the same time,
##  also the incidence relation for the orbits is computed.

(hash_list, edges) = GITFans.fan_traversal(fan_descr);


############################################################################
##  There are six maximal cones, up to G-symmetry.

length(hash_list)


############################################################################
##  We ask Polymake to create the incidence graph of the orbits,
##  and to visualize it.

intergraph = Polymake.graph.graph_from_edges(collect(edges));

Polymake.graph.visual(intergraph)


############################################################################
##  We translate the descriptions of the six maximal cones back
##  to cone objects ...

orbit_list = fan_descr[1];

result_cones = map(x -> Polymake.polytope.intersection(
                          GITFans.cones_from_bitlist(orbit_list, x)...), hash_list);


############################################################################
##  ... and expand their G-orbits

hom = GITFans.action_on_target(Q, G);
expanded = GITFans.orbit_cone_orbits(result_cones, hom);
orbit_lengths = map(length, expanded)


############################################################################
##  There are altogether 76 maximal cones.

sum(orbit_lengths)


############################################################################
##  The full intersection graph of the fan has (76 vertices and) 180 edges.
##  (A simpleminded visualization of this graph is not very enlightening.)

maxcones = vcat( expanded... );

full_edges = GITFans.edges_intersection_graph(maxcones, first(dims) - 1);

length(full_edges)

# full_intergraph = Polymake.graph.graph_from_edges(collect(full_edges));

# Polymake.graph.visual(full_intergraph)


############################################################################
##  We create the fan object in Polymake ...

rays_maxcones = [ [ convert( Vector{Rational{BigInt}}, cone.RAYS[i, :] )
                    for i in 1:size( cone.RAYS, 1 ) ]
                  for cone in maxcones ];

allrays = sort( collect( Set( vcat( rays_maxcones... ) ) ) )

index_maxcones = [ sort( [ findfirst( x -> x == v, allrays )-1 for v in rays ] )
                   for rays in rays_maxcones ]

inputrays = hcat( allrays...)'
fanobj = Polymake.fan.PolyhedralFan( INPUT_RAYS = inputrays, INPUT_CONES = index_maxcones )


############################################################################
##  ... and compute its F-vector.

fanobj.F_VECTOR


############################################################################
##  We check some statements from the paper.
##  There are 36 orbit cones of dimension 5, in 4 orbits.

oc = GITFans.orbit_cones( II, Q, G );
hom = GITFans.action_on_target( Q, G );
oco = GITFans.orbit_cone_orbits( oc, hom );
map( length, oco )

