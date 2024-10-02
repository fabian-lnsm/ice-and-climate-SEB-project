# Surface-energy balance in antarctica: Dependence of the surface climate on cloud coverage

## Project description

This project uses observations from a weather station in antarctica. Our goal is to predict the surface climate (i.e. surface temperature and melt)
at the given location in the case that the weather at the location is always cloudy. 

## Repository structure

The scripts are located in the src-folder where analysisSEB.ipynb is the main script to run.

The original data can be found in the data folder.

See the docs folder for the project instructions and the final report.

## Configuration

The code is written for Python 3.10 and makes use of the following modules: numpy 1.26, matplotlib 3.8, seaborn, sklearn and datetime.

The exact requirements to run the project are listed in the environment file which was created using conda.
The environment can be reproduced with the following command:

```
conda env create -f environment.yml
```



## Contribute

If you would like to suggest changes to the project, feel free to fork the project and make a pull request.

```
git fork https://github.com/fabian-lnsm/workshop-RC.git
```

## License

This project is licensed under the terms of the [MIT License](/LICENSE).
