# Silence R CMD check notes for the ggplot2 tidy-evaluation pronoun `.data`,
# which is only available at runtime (ggplot2 is an optional Suggests dependency).
utils::globalVariables(".data")
