#!/bin/bash
(cd solid && exec ./MolDyn --restart)
(cd gas && exec ./MolDyn --restart)
(cd liquid && exec ./MolDyn --restart)
