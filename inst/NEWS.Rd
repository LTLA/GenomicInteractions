\name{GenomicInteractionsNews}
\title{GenomicInteractions News}
\encoding{UTF-8}

\section{Version 2.0.0}{\itemize{
\item Completely new version of GenomicInteractions, with a redefined class based on the IndexedRelations package.
Users should updateObject() to update any serialized instances to this new definition.
\item Many high-level methods for manipulating GenomicInteractions instances have been moved to the fugi package.
For the time being, fugi is attached upon loading GenomicInteractions to preserve backwards compatibility.
This will not be the case in the next release so users should attach fugi explicitly. 
}}
