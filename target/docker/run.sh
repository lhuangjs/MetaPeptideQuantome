#!/bin/sh
#echo "Waiting for the eureka server to start on port $EUREKASERVER_PORT"
#while ! $(nc -z eureka-server $EUREKASERVER_PORT); do sleep 3; done
#echo "Eureka Server has started"
#
#echo "Waiting for the database server to start on port $DATABASESERVER_PORT ...................."
#while ! `nc -z database $DATABASESERVER_PORT`; echo $(nc -z database $DATABASESERVER_PORT 2>&1); do sleep 3; done
#echo "Database Server has started ......................."

echo "Starting metaproteomics-analysis-pipelines Service"
java -Djava.security.egd=file:/dev/./urandom \
  -Dserver.port=$SERVER_PORT \
  -Dspring.profiles.active=$PROFILE \
  -Dspring.datasource.url=jdbc:mysql://${DB_HOST}:${DB_PORT}/metaproteomics?useAffectedRows=true \
  -jar /usr/local/metaproteomics-analysis-pipelines/metaproteomics-analysis-pipelines-0.0.1-SNAPSHOT.jar
echo "Metaproteomics-analysis-pipelines Service has started"
