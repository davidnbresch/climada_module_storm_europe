climada_module_storm_europe (formerly named climada_module_ws_europe, changed 20151224)
===========================

European winter storm hazard event set

This repository contains an additional climada module, please install climada (the core module) first, see repository https://github.com/davidnbresch/climada

In order to grant core climada access to additional modules, create a folder ‘modules’ in the core climada folder and copy/move any additional modules into climada/modules, without 
'climada_module_' in the filename. You might also create a folder ‘parallel’ to climada (i.e. in the same folder as the core climada folder) and name it climada_modules to store additional modules, but this special setup is for developer use mainly).

E.g. if the addition module is named climada_module_MODULE_NAME, we should have
.../climada the core climada, with sub-folders as
.../climada/code
.../climada/data
.../climada/docs
and in there
.../climada/modules/MODULE_NAME with contents such as 
.../climada/modules/MODULE_NAME/code
.../climada/modules/MODULE_NAME/data
.../climada/modules/MODULE_NAME/docs
this way, climada sources all modules' code upon startup

see climada/docs/climada_manual.pdf to get started

copyright (c) 2014, David N. Bresch, david.bresch@gmail.com all rights reserved.
