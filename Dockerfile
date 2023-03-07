FROM continuumio/miniconda3

# set a directory for the app
WORKDIR /Users/CamillaTafuro/GitHub/microbiome/

# copy all the files to the container
COPY . .


# install dependencies
RUN pip install --no-cache-dir -r requirements.txt


# tell the port number the container should expose
EXPOSE 8501

# The code to run when container is started:
# run the command
CMD ["streamlit", "run", "./progetto_dashboard_microbioma.py"]
