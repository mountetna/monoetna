7c7
<       POSTGRES_DB: "polyphemus"
---
>       POSTGRES_DB: "${APP_NAME}"
15c15
<           - node.labels.volumes.polyphemus_staging_db == true
---
>           - node.labels.volumes.${APP_NAME}_${ENV}_db == true
24,25c24,25
<       - APP_NAME=polyphemus
<       - POLYPHEMUS_ENV=production
---
>       - APP_NAME=${APP_NAME}
>       - ${APP_ENV}=production
29a30,31
>       #- ETNA__PRODUCTION__HMAC_KEYS__METIS_FILE=/run/secrets/staging-metis-hmac
>       #- ETNA__PRODUCTION__HMAC_KEYS__MAGMA_FILE=/run/secrets/staging-metis-hmac
37c39,40
<     image: etnaagent/polyphemus:staging
---
> 
>     #user: ${APP_UID}:999 # metis:docker
43a47,49
>       #placement:
>         #constraints:
>           #- node.labels.volumes.${APP_NAME}_${ENV}_data == true
62c68
<           - polyphemus_app
---
>           - ${APP_NAME}_app
63a70,71
>       #- staging-metis-hmac
>       #- staging-magma-hmac
72,74c80,82
<       - APP_NAME=polyphemus
<       - POLYPHEMUS_ENV=production
<     image: etnaagent/polyphemus_app_fe:staging
---
>       - APP_NAME=${APP_NAME}
>       - ${APP_ENV}=production
>     image: etnaagent/${APP_NAME}_app_fe:${ENV}
79a88
>       #- ${DATA_SOURCE}/${APP_NAME}/:/app/data/${APP_NAME}:ro
88a98,99
>       restart_policy:
>         condition: none
91c102,103
<           - node.labels.volumes.big-share == true
---
>           - node.labels.volumes.${APP_NAME}_${ENV}_app_public == true
>           #- node.labels.volumes.${APP_NAME}_${ENV}_data == true
95,96c107,108
<       - APP_NAME=polyphemus
<       - POLYPHEMUS_ENV=production
---
>       - APP_NAME=${APP_NAME}
>       - ${APP_ENV}=production
99c111
<     image: etnaagent/polyphemus:staging
---
>     image: etnaagent/${APP_NAME}:${ENV}
112c124
<           - node.labels.volumes.big-share == true
---
>           - node.labels.volumes.${APP_NAME}_${ENV}_app_public == true
116,120d127
<     driver: local
<     driver_opts:
<       type: nfs
<       o: addr=mtetna-nfs1,nolock,noatime,suid,nodev,retry=16,proto=tcp,nfsvers=3
<       device: ":/cache/polyphemus/staging/app_public"
122c129
<     name: "polyphemus_staging_db"
---
>     name: "${APP_NAME}_${ENV}_db"
144a152,157
>   staging-metis-hmac:
>     external: true
>     name: ${STAGING_METIS_HMAC:-staging-metis-hmac}
>   staging-magma-hmac:
>     external: true
>     name: ${STAGING_MAGMA_HMAC:-staging-magma-hmac}
150c163,165
<     name: "polyphemus_staging_db_password"
---
>     name: "${APP_NAME}_${ENV}_db_password"
> 
>
