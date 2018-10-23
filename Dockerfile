FROM python:3-slim
RUN apt-get update
RUN apt-get install -yq gfortran
RUN apt-get install -yq liblapack-dev liblapacke-dev
RUN apt-get install -yq mongodb
ADD requirements.txt /app/
WORKDIR /app
RUN cat requirements.txt | xargs -n 1 -L 1 pip install
