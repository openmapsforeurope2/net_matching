#!/bin/bash
# gdb_conda_auto.sh
# Lance GDB avec Python embarqué fonctionnel, et restaure LD_LIBRARY_PATH Conda automatiquement à l'exécution de l'application

# Sauvegarde l'ancien LD_LIBRARY_PATH
export OLD_LD_LIBRARY_PATH=$LD_LIBRARY_PATH

# Forcer GDB à voir les libs système uniquement
export LD_LIBRARY_PATH=/usr/lib/x86_64-linux-gnu

# Lancer GDB avec un init file temporaire
TMP_GDBINIT=$(mktemp)

# Injecte un hook pour restaurer LD_LIBRARY_PATH à 'run'
cat << 'EOF' > $TMP_GDBINIT
define hook-run
    set env LD_LIBRARY_PATH=$OLD_LD_LIBRARY_PATH
end
EOF

# Lancer GDB en incluant le fichier init
gdb -x $TMP_GDBINIT "$@"

# Supprimer le fichier temporaire après fermeture
rm -f $TMP_GDBINIT