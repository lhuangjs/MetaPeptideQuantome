FROM openjdk:8-jdk-alpine
ENV ENV="/etc/profile"
RUN mkdir -p /data /usr/local/MetaPeptideQuantome /usr/local/MetaPeptideQuantome/software
ADD MetaPeptideQuantome-1.0-SNAPSHOT-spring-boot.jar /usr/local/MetaPeptideQuantome
ADD software /usr/local/MetaPeptideQuantome/software
RUN echo 'alias MPQ="java -jar /usr/local/MetaPeptideQuantome/MetaPeptideQuantome-1.0-SNAPSHOT-spring-boot.jar"' >> /etc/profile
