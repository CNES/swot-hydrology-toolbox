
WORKDIR=/work/ALT/swot/swotpub/SWOT_Simulator_data/SLCs_2_Products/po/20191106_bis_multitemp/
INDEX=0

debug=1

for tmp_simu in $(find $WORKDIR -maxdepth 6 -type f -name "pixcvec.nc"); do
    
    if [ ${tmp_simu:(-43):(-42)} -eq 3 ]; then
        if [ $debug -eq 1 ]; then
            cp ${tmp_simu} /work/ALT/swot/swotdev/desrochesd/floodplain/run/po/data/pixcvec/pixcvec${INDEX}.nc
        fi
    
        for tmp_simu in $(find ${tmp_simu:0:91} -maxdepth 6 -type f -name "pixel_cloud.nc"); do
            if [ ${tmp_simu:(-42):(-41)} -eq 3 ]; then
                if [ $debug -eq 1 ]; then
                    cp ${tmp_simu} /work/ALT/swot/swotdev/desrochesd/floodplain/run/po/data/pixc/pixel_cloud${INDEX}.nc
                fi
            fi
        done
        let INDEX=${INDEX}+1

    fi
done
