#This is an IRIT script and as such requires both math and irit import:
#
import math
import irit
#


# 
#  Curve Curve Intersection properties.
# 

view_mat2 = irit.sc( 1 )
irit.viewobj( view_mat2 )

# ######################################################

view_mat2 = irit.sc( 0.5 ) * irit.ty( 0.5 )
irit.viewobj( view_mat2 )

c1 = irit.cbspline( 2, irit.list( irit.ctlpt( irit.E3, (-2 ), 0, 0 ), \
                                  irit.ctlpt( irit.E2, 0, (-2 ) ), \
                                  irit.ctlpt( irit.E2, 2, 0 ) ), \
                    irit.list( irit.KV_OPEN ) )
c2 = irit.cbspline( 2, irit.list( irit.ctlpt( irit.E3, 2, (-2 ), 0 ), \
                                  irit.ctlpt( irit.E2, 0, 0 ), \
                                  irit.ctlpt( irit.E2, (-2 ), (-2 ) ) ), \
                    irit.list( irit.KV_OPEN ) )

ip = irit.ccintrprop( c1, c2, 1e-008 )
irit.printf( "number of intersection loops = %d\n", irit.list( irit.SizeOf( ip ) - 1 ) )

inter = irit.nth( ip, 2 )
irit.color( inter, 4 )
irit.adwidth( inter, 3 )

irit.printf( "area = %f\n", irit.list( irit.getattr( inter, "area" ) ) )

irit.interact( irit.list( c1, c2, inter ) )

# ################################

view_mat2 = irit.sc( 0.5 )
irit.viewobj( view_mat2 )

c1 = irit.cbspline( 4, irit.list( irit.ctlpt( irit.E3, 1, 0, 0 ), \
                                  irit.ctlpt( irit.E2, 1, 0.5523 ), \
                                  irit.ctlpt( irit.E2, 0.5523, 1 ), \
                                  irit.ctlpt( irit.E2, (-0.5523 ), 1 ), \
                                  irit.ctlpt( irit.E2, (-1 ), 0.5523 ), \
                                  irit.ctlpt( irit.E2, (-1 ), (-0.5523 ) ), \
                                  irit.ctlpt( irit.E2, (-0.5523 ), (-1 ) ), \
                                  irit.ctlpt( irit.E2, 0.5523, (-1 ) ), \
                                  irit.ctlpt( irit.E2, 1, (-0.5523 ) ), \
                                  irit.ctlpt( irit.E1, 1 ) ), \
                    irit.list( 0, 0, 0, 0, 1, 1, 2, 2, 3, 3, 4, 4, 4, 4 ) )
c2 = irit.cbezier( irit.list( irit.ctlpt( irit.E3, -1.4142, -1.4142, 0 ), \
                              irit.ctlpt( irit.E2, 1.4142, 1.4142 ) ) )

ip = irit.ccintrprop( c1, c2, 1e-008 )
irit.printf( "number of intersection loops = %d\n", irit.list( irit.SizeOf( ip ) - 1 ) )

inter = irit.nth( ip, 2 )
irit.color( inter, 4 )
irit.adwidth( inter, 3 )

irit.printf( "area = %f\n", irit.list( irit.getattr( inter, "area" ) ) )

all = irit.list( c1, c2, inter )

irit.interact( all )

irit.save( "cci1prop", all )

# ################################

view_mat2 = irit.sc( 0.2 ) * irit.ty( (-0.2 ) )
irit.viewobj( view_mat2 )

c1 = irit.cbspline( 3, irit.list( irit.ctlpt( irit.E3, 2.0018, 4.9791, 0 ), \
                                  irit.ctlpt( irit.E2, 4.452, 0.4139 ), \
                                  irit.ctlpt( irit.E2, (-5.4748 ), 0.4841 ), \
                                  irit.ctlpt( irit.E2, (-3.0089 ), 4.951 ) ), \
                    irit.list( irit.KV_OPEN ) )
c2 = irit.cbspline( 3, irit.list( irit.ctlpt( irit.E3, -3.0247, -1.9861, 0 ), \
                                  irit.ctlpt( irit.E2, (-5.4748 ), 2.5791 ), \
                                  irit.ctlpt( irit.E2, 4.452, 2.5089 ), \
                                  irit.ctlpt( irit.E2, 1.9861, (-1.958 ) ) ), \
                    irit.list( irit.KV_OPEN ) )

ip = irit.ccintrprop( c1, c2, 1e-008 )
irit.printf( "number of intersection loops = %d\n", irit.list( irit.SizeOf( ip ) - 1 ) )

inter = irit.nth( ip, 2 )
irit.color( inter, 4 )
irit.adwidth( inter, 3 )

irit.printf( "area = %f\n", irit.list( irit.getattr( inter, "area" ) ) )

all = irit.list( c1, c2, inter )

irit.interact( all )

irit.save( "cci2prop", all )

# ################################

view_mat2 = irit.sc( 0.2 ) * irit.ty( 0.3 )
irit.viewobj( view_mat2 )

c1 = irit.cbspline( 2, irit.list( irit.ctlpt( irit.E3, 4, (-3 ), 0 ), \
                                  irit.ctlpt( irit.E2, 2, (-1 ) ), \
                                  irit.ctlpt( irit.E2, 0, (-3 ) ), \
                                  irit.ctlpt( irit.E2, (-2 ), (-1 ) ), \
                                  irit.ctlpt( irit.E2, (-4 ), (-3 ) ) ), \
                    irit.list( irit.KV_OPEN ) )
c2 = irit.cbspline( 2, irit.list( irit.ctlpt( irit.E3, (-4 ), (-1 ), 0 ), \
                                  irit.ctlpt( irit.E2, (-2 ), (-3 ) ), \
                                  irit.ctlpt( irit.E2, 0, (-1 ) ), \
                                  irit.ctlpt( irit.E2, 2, (-3 ) ), \
                                  irit.ctlpt( irit.E2, 4, (-1 ) ) ), \
                    irit.list( irit.KV_OPEN ) )

ip = irit.ccintrprop( c1, c2, 1e-008 )
irit.printf( "number of intersection loops = %d\n", irit.list( irit.SizeOf( ip ) - 1 ) )

inter = irit.list( irit.nth( ip, 2 ), irit.nth( ip, 3 ) )
irit.color( inter, 4 )
irit.adwidth( inter, 3 )

irit.printf( "area = (%f, %f)\n", irit.list( irit.getattr( irit.nth( ip, 2 ), "area" ), irit.getattr( irit.nth( ip, 3 ), "area" ) ) )

all = irit.list( c1, c2, inter )

irit.interact( all )

irit.save( "cci3prop", all )

# ################################

view_mat2 = irit.sc( 0.3 ) * irit.ty( (-0.2 ) )
irit.viewobj( view_mat2 )

c1 = irit.cbspline( 4, irit.list( irit.ctlpt( irit.E3, (-2.921 ), (-0.2957 ), 0 ), \
                                  irit.ctlpt( irit.E2, (-2.4402 ), (-1.2055 ) ), \
                                  irit.ctlpt( irit.E2, (-1.4524 ), (-0.6661 ) ), \
                                  irit.ctlpt( irit.E2, (-0.991 ), 0.2697 ), \
                                  irit.ctlpt( irit.E2, (-0.4841 ), 0.3282 ), \
                                  irit.ctlpt( irit.E2, 0.5361, (-0.3477 ) ), \
                                  irit.ctlpt( irit.E2, 1.4979, (-1.0105 ) ), \
                                  irit.ctlpt( irit.E2, 2.2127, (-0.874 ) ), \
                                  irit.ctlpt( irit.E2, 2.5636, 0.0227 ), \
                                  irit.ctlpt( irit.E2, 1.7058, 1.199 ), \
                                  irit.ctlpt( irit.E2, 0.6661, 1.4329 ), \
                                  irit.ctlpt( irit.E2, 0.8091, 1.7383 ), \
                                  irit.ctlpt( irit.E2, 1.4784, 1.9138 ), \
                                  irit.ctlpt( irit.E2, 2.1477, 1.5824 ), \
                                  irit.ctlpt( irit.E2, 2.7651, 1.8358 ) ), \
                    irit.list( irit.KV_OPEN ) )
c2 = irit.cbspline( 4, irit.list( irit.ctlpt( irit.E3, 2.6806, 1.4264, 0 ), \
                                  irit.ctlpt( irit.E2, 2.1672, 1.9463 ), \
                                  irit.ctlpt( irit.E2, 2.1022, 1.4654 ), \
                                  irit.ctlpt( irit.E2, 0.5946, 1.8618 ), \
                                  irit.ctlpt( irit.E2, 0.3022, 1.2964 ), \
                                  irit.ctlpt( irit.E2, 0.7896, 0.8415 ), \
                                  irit.ctlpt( irit.E2, 2.3687, 0.8415 ), \
                                  irit.ctlpt( irit.E2, 2.8106, 0.3022 ), \
                                  irit.ctlpt( irit.E2, 2.6871, (-0.4191 ) ), \
                                  irit.ctlpt( irit.E2, 1.9723, (-0.2892 ) ), \
                                  irit.ctlpt( irit.E2, 1.2314, (-0.1722 ) ), \
                                  irit.ctlpt( irit.E2, 0.5621, (-0.5426 ) ), \
                                  irit.ctlpt( irit.E2, (-0.0552 ), (-1.1145 ) ), \
                                  irit.ctlpt( irit.E2, (-0.5881 ), (-0.7831 ) ), \
                                  irit.ctlpt( irit.E2, (-1.4979 ), 0.0032 ), \
                                  irit.ctlpt( irit.E2, (-2.1087 ), 0.0162 ), \
                                  irit.ctlpt( irit.E2, (-2.8886 ), (-0.7766 ) ) ), \
                    irit.list( irit.KV_OPEN ) )

ip = irit.ccintrprop( c1, c2, 1e-008 )
irit.printf( "number of intersection loops = %d\n", irit.list( irit.SizeOf( ip ) - 1 ) )

inter = irit.list( irit.nth( ip, 2 ), irit.nth( ip, 3 ), irit.nth( ip, 4 ), irit.nth( ip, 5 ) )
irit.color( inter, 4 )
irit.adwidth( inter, 3 )

all = irit.list( c1, c2, inter )

irit.printf( "area = (%f, %f, %f, %f)\n", irit.list( irit.getattr( irit.nth( ip, 2 ), "area" ), irit.getattr( irit.nth( ip, 3 ), "area" ), irit.getattr( irit.nth( ip, 4 ), "area" ), irit.getattr( irit.nth( ip, 5 ), "area" ) ) )

dadpt = irit.nth( ip, 1 )
n = irit.SizeOf( dadpt )

irit.printf( "d area / d control points of 1st curve:\n", irit.nil(  ) )
i = 1
while ( i <= irit.SizeOf( c1 ) ):
    irit.printf( "\tpt(%d)x = %f  pt(%d)y = %f\n", irit.list( i, irit.nth( dadpt, i * 2 - 1 ), i, irit.nth( dadpt, i * 2 ) ) )
    i = i + 1

irit.printf( "d area / d control points of 2nd curve:\n", irit.nil(  ) )
j = 1
while ( j <= irit.SizeOf( c2 ) ):
    irit.printf( "\tpt(%d)x = %f  pt(%d)y = %f\n", irit.list( j, irit.nth( dadpt, ( i + j - 1 ) * 2 - 1 ), j, irit.nth( dadpt, ( i + j - 1 ) * 2 ) ) )
    j = j + 1

irit.interact( all )

irit.save( "cci4prop", all )

# ################################
#  Local intersections - 
#      some dAdPt should be zero!
# ################################

view_mat2 = irit.sc( 0.25 ) * irit.tx( (-0.2 ) ) * irit.ty( (-0.2 ) )
irit.viewobj( view_mat2 )

c1 = irit.cbspline( 4, irit.list( irit.ctlpt( irit.E3, (-2.921 ), (-0.2957 ), 0 ), \
                                  irit.ctlpt( irit.E2, (-2.4402 ), (-1.2055 ) ), \
                                  irit.ctlpt( irit.E2, (-1.4524 ), (-0.6661 ) ), \
                                  irit.ctlpt( irit.E2, (-0.991 ), 0.2697 ), \
                                  irit.ctlpt( irit.E2, (-0.4841 ), 0.3282 ), \
                                  irit.ctlpt( irit.E2, 0.5361, (-0.3477 ) ), \
                                  irit.ctlpt( irit.E2, 1.4979, (-1.0105 ) ), \
                                  irit.ctlpt( irit.E2, 2.2127, (-0.874 ) ), \
                                  irit.ctlpt( irit.E2, 2.5636, 0.0227 ), \
                                  irit.ctlpt( irit.E2, 1.7058, 1.199 ), \
                                  irit.ctlpt( irit.E2, 0.6661, 1.4329 ), \
                                  irit.ctlpt( irit.E2, 0.8091, 1.7383 ), \
                                  irit.ctlpt( irit.E2, 1.4784, 1.9138 ), \
                                  irit.ctlpt( irit.E2, 2.1477, 1.5824 ), \
                                  irit.ctlpt( irit.E2, 2.7651, 1.8358 ) ), \
                    irit.list( irit.KV_OPEN ) )
c2 = irit.cbspline( 4, irit.list( irit.ctlpt( irit.E3, 2.5882, 3.5517, 0 ), \
                                  irit.ctlpt( irit.E2, 1.069, 3.4251 ), \
                                  irit.ctlpt( irit.E2, 1.5613, 2.5108 ), \
                                  irit.ctlpt( irit.E2, 4.0088, 2.7499 ), \
                                  irit.ctlpt( irit.E2, 4.5152, 2.7077 ), \
                                  irit.ctlpt( irit.E2, 3.7697, 2.0607 ), \
                                  irit.ctlpt( irit.E2, 3.9666, 1.3996 ), \
                                  irit.ctlpt( irit.E2, 3.5587, 0.4853 ), \
                                  irit.ctlpt( irit.E2, 3.3618, (-1.0057 ) ), \
                                  irit.ctlpt( irit.E2, 0.9284, (-2.4686 ) ), \
                                  irit.ctlpt( irit.E2, (-0.3798 ), (-2.2154 ) ), \
                                  irit.ctlpt( irit.E2, (-0.647 ), (-1.5684 ) ), \
                                  irit.ctlpt( irit.E2, (-0.0552 ), (-1.1145 ) ), \
                                  irit.ctlpt( irit.E2, (-0.5881 ), (-0.7831 ) ), \
                                  irit.ctlpt( irit.E2, (-1.4979 ), 0.0032 ), \
                                  irit.ctlpt( irit.E2, (-2.1087 ), 0.0162 ), \
                                  irit.ctlpt( irit.E2, (-2.8886 ), (-0.7766 ) ) ), \
                    irit.list( irit.KV_OPEN ) )

ip = irit.ccintrprop( c1, c2, 1e-008 )
irit.printf( "number of intersection loops = %d\n", irit.list( irit.SizeOf( ip ) - 1 ) )

inter = irit.list( irit.nth( ip, 2 ) )
irit.color( inter, 4 )
irit.adwidth( inter, 3 )

all = irit.list( c1, c2, inter )

irit.printf( "area = (%f)\n", irit.list( irit.getattr( irit.nth( ip, 2 ), "area" ) ) )

dadpt = irit.nth( ip, 1 )
n = irit.SizeOf( dadpt )

irit.printf( "d area / d control points of 1st curve:\n", irit.nil(  ) )
i = 1
while ( i <= irit.SizeOf( c1 ) ):
    irit.printf( "\tpt(%d)x = %f  pt(%d)y = %f\n", irit.list( i, irit.nth( dadpt, i * 2 - 1 ), i, irit.nth( dadpt, i * 2 ) ) )
    i = i + 1

irit.printf( "d area / d control points of 2nd curve:\n", irit.nil(  ) )
j = 1
while ( j <= irit.SizeOf( c2 ) ):
    irit.printf( "\tpt(%d)x = %f  pt(%d)y = %f\n", irit.list( j, irit.nth( dadpt, ( i + j - 1 ) * 2 - 1 ), j, irit.nth( dadpt, ( i + j - 1 ) * 2 ) ) )
    j = j + 1

irit.interact( all )

irit.save( "cci5prop", all )

# ################################

irit.free( dadpt )
irit.free( c1 )
irit.free( c2 )
irit.free( ip )
irit.free( inter )
irit.free( all )
irit.free( view_mat2 )

