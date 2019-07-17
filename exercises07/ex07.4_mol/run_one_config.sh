#!/bin/bash
(cd solid && exec ./MolDyn --one-config)
(cd gas && exec ./MolDyn --one-config)
(cd liquid && exec ./MolDyn --one-config)
