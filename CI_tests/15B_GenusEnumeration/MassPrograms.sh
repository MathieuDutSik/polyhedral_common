# The namelists of GENUS_IsometryMass and GENUS_ComponentMass, written from
# key=value arguments into a temporary file, the program run on it and its
# result printed on stdout:
#
#   isometry_mass GramFile=isom_E8.txt GeneratorsFile=isom_g.txt [mode=check]
#   component_mass field=abelian_real output=local_factors n=5 H=4 ...
#
# The keys go to their block (DATA, METHOD, SYSTEM) by name, the integer
# ones unquoted; OutFile is stdout unless given. Sourced by the test scripts.

MASS_SRC=../../src_genera

run_mass_program() {
    prog=$1
    shift
    data=""
    method=""
    system=' OutFile = "stdout"\n'
    for kv in "$@"; do
        k=${kv%%=*}
        v=${kv#*=}
        case $k in
            n|rank|hecke_X|hecke_prec_bits) val=$v ;;
            *) val="\"$v\"" ;;
        esac
        case $k in
            mode|check_kind|lvalue|output|hecke_X|hecke_prec_bits)
                method="$method $k = $val\n" ;;
            OutFile) system=" OutFile = $val\n" ;;
            *) data="$data $k = $val\n" ;;
        esac
    done
    nml=$(mktemp)
    {
        printf '&DATA\n%b/\n\n&METHOD\n%b/\n' "$data" "$method"
        if [ "$prog" = GENUS_IsometryMass ]; then
            printf '\n&SYSTEM\n%b/\n' "$system"
        fi
    } > "$nml"
    status=0
    $MASS_SRC/$prog "$nml" || status=$?
    rm -f "$nml"
    return $status
}

isometry_mass() {
    run_mass_program GENUS_IsometryMass "$@"
}

component_mass() {
    run_mass_program GENUS_ComponentMass "$@"
}
