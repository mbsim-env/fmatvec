// load mathjax-node module or return with a warning message
try {
  var mj=require("mathjax-node");
}
catch(ex) {
  console.log("Warning: Skipping MathJax check. nodejs module mathjax-node not found in NODE_PATH.");
  process.exitCode=0;
  return;
}

console.log("MathJax checking using nodejs module mathjax-node of files:");
for(var i in process.argv) {
  if(i==0 || i==1)
    continue;
  console.log(process.argv[i]);
}

// load standard modules
var fs=require('fs');

// configuration

var inlineMath=[['\\(','\\)']];
var displayMath=[['\\[','\\]']];

mj.config({
  displayErrors: false,
  undefinedCharError: true,
  MathJax: { 
    extensions: ["tex2jax.js", "TeX/AMSmath.js", "TeX/AMSsymbols.js"],
    jax: ["input/TeX","output/HTML-CSS"],
    tex2jax: {
      inlineMath: inlineMath,
      displayMath: displayMath,
    },
    TeX: {
      Macros: {
        dd: "\\mathrm{d}",
        vM: "\\boldsymbol{M}",
        vK: "\\boldsymbol{K}",
        vD: "\\boldsymbol{D}",
        vW: "\\boldsymbol{W}",
        vV: "\\boldsymbol{V}",
        vJ: "\\boldsymbol{J}",
        vF: "\\boldsymbol{F}",
        vq: "\\boldsymbol{q}",
        vu: "\\boldsymbol{u}",
        vx: "\\boldsymbol{x}",
        vg: "\\boldsymbol{g}",
        vz: "\\boldsymbol{z}",
        vh: "\\boldsymbol{h}",
        vs: "\\boldsymbol{s}",
        vr: "\\boldsymbol{r}",
        vv: "\\boldsymbol{v}",
        vLambda: "\\boldsymbol{\\Lambda}",
        vlambda: "\\boldsymbol{\\lambda}",
      },
    },
  }
});

// init mathjax
mj.start();

// global mathjax state (used as a cache to improve performance)
var mjState={};

// loop over all arguments
for(var i in process.argv) {
  if(i==0 || i==1)
    continue;
  processFileOrDir(process.argv[i]);
}

// handle file or directory path
function processFileOrDir(path) {
  var s=fs.lstatSync(path);
  if(s.isFile())
    processFile(path);
  if(s.isDirectory())
    processPath(path);
}

// recursively walk directories
function processPath(path) {
  var files=fs.readdirSync(path);
  for(var i in files)
    processFileOrDir(path+"/"+files[i]);
}

// process a file path
function processFile(path) {
  // skip none .html files
  if(!path.endsWith('.html'))
    return;

  // read file content
  var data=fs.readFileSync(path, 'utf8');
  // check all inline math deliminiter
  for(var i in inlineMath)
    checkEqn(path, data, inlineMath[i][0], inlineMath[i][1]);
  // check all display math deliminiter
  for(var i in displayMath)
    checkEqn(path, data, displayMath[i][0], displayMath[i][1]);
}

// check math given the the deliminiter
function checkEqn(filename, data, startTag, endTag) {
  var s;
  // loop over all geliminiters
  while((s=data.indexOf(startTag))>=0) {
    data=data.substr(s+startTag.length);
    var e=data.indexOf(endTag);
    if(e<0) {
      console.log(filename+": Opening MathJax tag found but no MathJax end tag found.");
      process.exitCode=1;
      return;
    }
    // get the equation between the deliminters
    var eqnStr=data.substr(0, e);
    // predefined entities
    eqnStr=eqnStr.replace(/&quot;/g, '"');
    eqnStr=eqnStr.replace(/&amp;/g , '&');
    eqnStr=eqnStr.replace(/&apos;/g, "'");
    eqnStr=eqnStr.replace(/&lt;/g  , '<');
    eqnStr=eqnStr.replace(/&gt;/g  , '>');
    // numeric character reference 
    // MISSING: replace &#nnn; with String.fromCodePoint(parseInt(nnn))
    // MISSING: replace &#xhhh; with String.fromCodePoint(parseInt("0x"+hex))

    // process the rest of the file content
    data=data.substr(e+endTag.length);

    // process the equation with mathjax
    mj.typeset({
      math: eqnStr,
      format: "TeX",
      state: mjState,
    }, function(formula, data) {
      if(data.errors) {
        // on errors print the mathjax error message
        for(var i in data.errors) {
          process.exitCode=1;
          console.log(filename+": "+formula+": "+data.errors[i]);
        }
        mjState=data.state;
      }
    }.bind(null, eqnStr));
  }
}
