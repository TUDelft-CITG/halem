# Start with pyramid app image
FROM continuumio/miniconda3

# Install conda stuff first
RUN conda install numpy pandas nomkl pyproj
# Then install rest via pip
RUN pip install pint Flask

ADD . /Route_optimization_in_dynamic_currents
WORKDIR /Route_optimization_in_dynamic_currents

# Install the application
RUN pip install -e .

# expose port 5000
EXPOSE 5000
# Serve on port 5000
CMD halem serve --port 5000
