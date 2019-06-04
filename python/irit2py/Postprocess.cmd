Rem
Rem post process an irit2py result - replacing all left overs...
Rem

sed -e "s/irit.kv_open/irit.KV_OPEN/" < %1 | sed -e "s/irit.sizeof/irit.SizeOf/"
