window.MathJax = {
  tex2jax: {
    inlineMath: [ ["\\(","\\)"] ],
    displayMath: [ ["\\[","\\]"] ]
  },
  TeX: {
    TagSide: "right",
    TagIndent: ".8em",
    MultLineWidth: "80%",
    equationNumbers: {
      autoNumber: "AMS",
    },
    unicode: {
      fonts: "STIXGeneral,'Arial Unicode MS'"
    },
    macros: {
      RR: "{\\bf R}",
      bold: ["{\\bf #1}", 1]
    }
  },
  displayAlign: "center",
  showProcessingMessages: false,
  messageStyle: "none"
};
