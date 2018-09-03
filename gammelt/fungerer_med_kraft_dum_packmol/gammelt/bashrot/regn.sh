function py {
    python -c "from math import *; x=$1; x=int(x) if int(x)==x else x; print(str(x))"
}
