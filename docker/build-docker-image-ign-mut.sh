
#!/bin/sh
if [ -d "../config_tmp" ]; then
    echo "WARNING : folder config_tmp deletion"
    mv ../config_tmp/db_conf* ../config
    rm -r ../config_tmp
fi

mkdir ../config_tmp
mv ../config/db_conf* ../config_tmp

./build-docker-image.sh

mv ../config_tmp/db_conf* ../config
rm -rf ../config_tmp