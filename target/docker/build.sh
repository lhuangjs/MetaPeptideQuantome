mvn clean package -Dmaven.test.skip=true
cd target/docker
cp ../MetaPeptideQuantome-1.0-SNAPSHOT-spring-boot.jar ./
docker build -t huangjs2017/meta-peptide-quantome:1.0 .