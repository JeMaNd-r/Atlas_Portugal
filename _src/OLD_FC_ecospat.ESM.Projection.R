ecospat.ESM.Projection.fixed <- function (ESM.modeling.output, new.env, name.env = NULL, parallel = FALSE, 
          cleanup = FALSE) 
{
  iniwd <- getwd()
  on.exit(setwd(iniwd))
  setwd(ESM.modeling.output$wd)
  models <- ESM.modeling.output$models
  models. <- ESM.modeling.output$models.
  mymodels <- ESM.modeling.output$mymodels
  data <- ESM.modeling.output$data
  combinations <- combn(colnames(ESM.modeling.output$data@data.env.var), 
                        2)
  which.biva <- ESM.modeling.output$which.biva
  NbRunEval <- ESM.modeling.output$NbRunEval
  modeling.id <- ESM.modeling.output$modeling.id
  if (is.matrix(new.env)) {
    new.env <- as.data.frame(new.env)
  }
  if (is.null(name.env)) 
    name.env <- deparse(substitute(new.env))
  if (parallel == FALSE) {
    for (k in 1:length(mymodels)) {
      mymodel <- mymodels[[k]]
      dir.create(path = paste0("./", "ESM.BIOMOD.", k, 
                               "/proj_", paste(name.env, "ESM.BIOMOD", k, modeling.id, 
                                               sep = ".")))
      if (is.character(mymodel)) {
        (next)()
      }
      if (is.data.frame(new.env)) {
        newdata = new.env[, colnames(new.env) %in% combinations[, 
                                                                k]]
        if ("PA.table" %in% slotNames(data)) {
          models.chosen = grep("allRun", mymodel@models.computed[-grep("allData", 
                                                                       mymodel@models.computed)], value = TRUE)
        }
        else {
          models.chosen = grep("allRun", mymodel@models.computed, 
                               value = TRUE)
        }
        for (i in 1:length(models)) {
          modelToProject <- get(biomod2::BIOMOD_LoadModels(mymodel, 
                                                           full.name = grep(paste0("_", models[i]), 
                                                                            models.chosen, value = TRUE)))
          if (models[i] == "MAXENT") {
            map <- biomod2::predict(modelToProject, 
                                    newdata = newdata, temp_workdir = modelToProject@model_output_dir, 
                                    overwrite = T, on_0_1000 = TRUE, omit.na = TRUE)
          }
          else {
            map <- biomod2::predict(modelToProject, 
                                    newdata = newdata, overwrite = T, on_0_1000 = TRUE, 
                                    omit.na = TRUE)
          }
          save(map, file = file.path(paste0("./", "ESM.BIOMOD.", 
                                            k, "/proj_", paste(name.env, "ESM.BIOMOD", 
                                                               k, modeling.id, sep = "."), "/proj_", 
                                            paste(name.env, models[i], "ESM.BIOMOD", 
                                                  k, modeling.id, sep = "."), ".RData")))
        }
      }
      if (inherits(new.env, "RasterStack")) {
        new.env = terra::rast(new.env)
      }
      if (inherits(new.env, "SpatRaster")) {
        newdata = newdata = subset(new.env, combinations[, 
                                                         k])
        if ("PA.table" %in% slotNames(data)) {
          models.chosen = grep("allRun", mymodel@models.computed[-grep("allData", 
                                                                       mymodel@models.computed)], value = TRUE)
        }
        else {
          models.chosen = grep("allRun", mymodel@models.computed, 
                               value = TRUE)
        }
        for (i in 1:length(models)) {
          modelToProject <- biomod2::BIOMOD_LoadModels(mymodel, 
                                                           full.name = grep(paste0("_", models[i]), 
                                                                            models.chosen, value = TRUE))
          if (models[i] == "MAXENT") {
            map <- biomod2::predict(modelToProject, 
                                    newdata = newdata, temp_workdir = modelToProject@model_output_dir, 
                                    overwrite = T, on_0_1000 = TRUE, omit.na = TRUE)
          }
          else {
            map <- biomod2::predict(modelToProject, 
                                    newdata = newdata, overwrite = T, on_0_1000 = TRUE, 
                                    omit.na = TRUE)
          }
          terra::writeRaster(map, paste0("./", "ESM.BIOMOD.", 
                                         k, "/proj_", paste(name.env, "ESM.BIOMOD", 
                                                            k, modeling.id, sep = "."), "/proj_", 
                                         paste(name.env, models[i], "ESM.BIOMOD", 
                                               k, modeling.id, sep = "."), ".tif"), overwrite = TRUE)
        }
      }
    }
  }
  if (parallel == TRUE) {
    if (inherits(new.env, "RasterStack")) {
      new.env = terra::rast(new.env)
    }
    if (inherits(new.env, "SpatRaster")) {
      new.env <- terra::wrap(new.env)
    }
    for (g in 1:length(mymodels)) {
      dir.create(path = paste0("./", "ESM.BIOMOD.", g, 
                               "/proj_", paste(name.env, "ESM.BIOMOD", g, modeling.id, 
                                               sep = ".")))
    }
    foreach(k = 1:length(mymodels), .packages = c("biomod2", 
                                                  "raster", "terra", "base")) %dopar% {
                                                    mymodel <- mymodels[[k]]
                                                    if (!(is.character(mymodel))) {
                                                      i = 0
                                                      setwd(ESM.modeling.output$wd)
                                                      if (is.data.frame(new.env)) {
                                                        newdata = new.env[, colnames(new.env) %in% 
                                                                            combinations[, k]]
                                                        if ("PA.table" %in% slotNames(data)) {
                                                          models.chosen = grep("allRun", mymodel@models.computed[-grep("allData", 
                                                                                                                       mymodel@models.computed)], value = TRUE)
                                                        }
                                                        else {
                                                          models.chosen = grep("allRun", mymodel@models.computed, 
                                                                               value = TRUE)
                                                        }
                                                        for (i in 1:length(models)) {
                                                          modelToProject <- biomod2::BIOMOD_LoadModels(mymodel, 
                                                                                                           full.name = grep(paste0("_", models[i]), 
                                                                                                                            models.chosen, value = TRUE))
                                                          if (models[i] == "MAXENT") {
                                                            map <- biomod2::predict(modelToProject, 
                                                                                    newdata = newdata, temp_workdir = modelToProject@model_output_dir, 
                                                                                    overwrite = T, on_0_1000 = TRUE, omit.na = TRUE)
                                                          }
                                                          else {
                                                            map <- biomod2::predict(modelToProject, 
                                                                                    newdata = newdata, overwrite = T, on_0_1000 = TRUE, 
                                                                                    omit.na = TRUE)
                                                          }
                                                          save(map, file = file.path(paste0("./", 
                                                                                            "ESM.BIOMOD.", k, "/proj_", paste(name.env, 
                                                                                                                              "ESM.BIOMOD", k, modeling.id, sep = "."), 
                                                                                            "/proj_", paste(name.env, models[i], "ESM.BIOMOD", 
                                                                                                            k, modeling.id, sep = "."), ".RData")))
                                                        }
                                                      }
                                                      else {
                                                        new.env <- terra::rast(new.env)
                                                        newdata <- terra::subset(new.env, combinations[, 
                                                                                                       k])
                                                        new.env <- terra::wrap(new.env)
                                                        if ("PA.table" %in% slotNames(data)) {
                                                          models.chosen = grep("allRun", mymodel@models.computed[-grep("allData", 
                                                                                                                       mymodel@models.computed)], value = TRUE)
                                                        }
                                                        else {
                                                          models.chosen = grep("allRun", mymodel@models.computed, 
                                                                               value = TRUE)
                                                        }
                                                        for (i in 1:length(models)) {
                                                          modelToProject <- biomod2::BIOMOD_LoadModels(mymodel, 
                                                                                                           full.name = grep(paste0("_", models[i]), 
                                                                                                                            models.chosen, value = TRUE))
                                                          if (models[i] == "MAXENT") {
                                                            map <- biomod2::predict(modelToProject, 
                                                                                    newdata = newdata, temp_workdir = modelToProject@model_output_dir, 
                                                                                    on_0_1000 = TRUE, omit.na = TRUE, overwrite = T)
                                                          }
                                                          else {
                                                            map <- biomod2::predict(modelToProject, 
                                                                                    newdata = newdata, on_0_1000 = TRUE, 
                                                                                    omit.na = TRUE, overwrite = T)
                                                          }
                                                          terra::writeRaster(map, file = file.path(paste0(ESM.modeling.output$wd, 
                                                                                                          "/", "ESM.BIOMOD.", k, "/proj_", paste(name.env, 
                                                                                                                                                 "ESM.BIOMOD", k, modeling.id, sep = "."), 
                                                                                                          "/proj_", paste(name.env, models[i], "ESM.BIOMOD", 
                                                                                                                          k, modeling.id, sep = "."), ".tif")), 
                                                                             overwrite = T)
                                                        }
                                                      }
                                                    }
                                                  }
    if (!is.data.frame(new.env)) {
      new.env <- terra::rast(new.env)
    }
  }
  if (cleanup != FALSE) {
    removeTmpFiles(h = cleanup)
  }
  output <- list(proj.name = name.env, modeling.id = modeling.id, 
                 models. = grep(modeling.id, gtools::mixedsort(list.files(getwd(), 
                                                                          "models.out", recursive = TRUE, full.names = TRUE)), 
                                value = TRUE), models = models, pred.biva = grep(modeling.id, 
                                                                                 gtools::mixedsort(list.files(getwd(), paste("proj_", 
                                                                                                                             name.env, sep = ""), recursive = TRUE, full.names = TRUE)), 
                                                                                 value = TRUE), NbRunEval = NbRunEval, name.env = name.env, 
                 new.env.raster = inherits(new.env, "SpatRaster"), wd = getwd(), 
                 which.biva = which.biva)
  save(output, file = paste("ESM_Projections", name.env, modeling.id, 
                            "out", sep = "."))
  setwd(iniwd)
  return(output)
}
