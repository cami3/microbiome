# Define base image
FROM continuumio/miniconda3
 
# Set working directory for the project
WORKDIR /app

LABEL maintainer="Camilla Tafuro"

# copy or clone all the files to the container
RUN git clone https://github.com/cami3/microbiome.git .

# Create Conda environment from the YAML file
COPY qiime2-2023.2-py38-linux-conda.yml .
RUN conda env create --name my_env -f qiime2-2023.2-py38-linux-conda.yml
 

# Override default shell and use bash
SHELL ["conda", "run", "-n", "my_env", "/bin/bash", "-c"]


# Activate Conda environment and check if it is working properly
RUN echo "Making sure qiime2 is installed correctly..."
RUN python -c "import qiime2"

# Install pip requirements
ADD ./requirements.txt .
RUN pip3 install -r ./requirements.txt

# tell the port number the container should expose
EXPOSE 8501

HEALTHCHECK CMD curl --fail http://localhost:8501/_stcore/health

# Python program to run in the container

# The code to run when container is started:
# run the command
ENTRYPOINT ["conda", "run", "-n", "my_env", "streamlit", "run", "progetto_dashboard_microbioma.py", "--server.port=8501", "--server.address=0.0.0.0"]

