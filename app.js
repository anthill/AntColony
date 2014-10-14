(function e(t,n,r){function s(o,u){if(!n[o]){if(!t[o]){var a=typeof require=="function"&&require;if(!u&&a)return a(o,!0);if(i)return i(o,!0);var f=new Error("Cannot find module '"+o+"'");throw f.code="MODULE_NOT_FOUND",f}var l=n[o]={exports:{}};t[o][0].call(l.exports,function(e){var n=t[o][1][e];return s(n?n:e)},l,l.exports,e,t,n,r)}return n[o].exports}var i=typeof require=="function"&&require;for(var o=0;o<r.length;o++)s(r[o]);return s})({"/Users/vallette/ANTS/AntColony/index.js":[function(require,module,exports){
'use strict';

var initRendering = require('./src/rendering.js');

module.exports = function(containerElement){
    initRendering(containerElement);
    var points = require('./src/initializePoints.js');
    var edges = require('./src/createEdges.js');
};
},{"./src/createEdges.js":"/Users/vallette/ANTS/AntColony/src/createEdges.js","./src/initializePoints.js":"/Users/vallette/ANTS/AntColony/src/initializePoints.js","./src/rendering.js":"/Users/vallette/ANTS/AntColony/src/rendering.js"}],"/Users/vallette/ANTS/AntColony/node_modules/delaunay-triangulate/node_modules/incremental-convex-hull/ich.js":[function(require,module,exports){
"use strict"

//High level idea:
// 1. Use Clarkson's incremental construction to find convex hull
// 2. Point location in triangulation by jump and walk

module.exports = incrementalConvexHull

var orient = require("robust-orientation")
var compareCell = require("simplicial-complex").compareCells

function compareInt(a, b) {
  return a - b
}

function Simplex(vertices, adjacent, boundary) {
  this.vertices = vertices
  this.adjacent = adjacent
  this.boundary = boundary
  this.lastVisited = -1
}

Simplex.prototype.flip = function() {
  var t = this.vertices[0]
  this.vertices[0] = this.vertices[1]
  this.vertices[1] = t
  var u = this.adjacent[0]
  this.adjacent[0] = this.adjacent[1]
  this.adjacent[1] = u
}

function GlueFacet(vertices, cell, index) {
  this.vertices = vertices
  this.cell = cell
  this.index = index
}

function compareGlue(a, b) {
  return compareCell(a.vertices, b.vertices)
}

function bakeOrient(d) {
  var code = ["function orient(){var tuple=this.tuple;return test("]
  for(var i=0; i<=d; ++i) {
    if(i > 0) {
      code.push(",")
    }
    code.push("tuple[", i, "]")
  }
  code.push(")}return orient")
  var proc = new Function("test", code.join(""))
  var test = orient[d+1]
  if(!test) {
    test = orient
  }
  return proc(test)
}

var BAKED = []

function Triangulation(dimension, vertices, simplices) {
  this.dimension = dimension
  this.vertices = vertices
  this.simplices = simplices
  this.interior = simplices.filter(function(c) {
    return !c.boundary
  })

  this.tuple = new Array(dimension+1)
  for(var i=0; i<=dimension; ++i) {
    this.tuple[i] = this.vertices[i]
  }

  var o = BAKED[dimension]
  if(!o) {
    o = BAKED[dimension] = bakeOrient(dimension)
  }
  this.orient = o
}

var proto = Triangulation.prototype

//Degenerate situation where we are on boundary, but coplanar to face
proto.handleBoundaryDegeneracy = function(cell, point) {
  var d = this.dimension
  var n = this.vertices.length - 1
  var tuple = this.tuple
  var verts = this.vertices

  //Dumb solution: Just do dfs from boundary cell until we find any peak, or terminate
  var toVisit = [ cell ]
  cell.lastVisited = -n
  while(toVisit.length > 0) {
    cell = toVisit.pop()
    var cellVerts = cell.vertices
    var cellAdj = cell.adjacent
    for(var i=0; i<=d; ++i) {
      var neighbor = cellAdj[i]
      if(!neighbor.boundary || neighbor.lastVisited <= -n) {
        continue
      }
      var nv = neighbor.vertices
      for(var j=0; j<=d; ++j) {
        var vv = nv[j]
        if(vv < 0) {
          tuple[j] = point
        } else {
          tuple[j] = verts[vv]
        }
      }
      var o = this.orient()
      if(o > 0) {
        return neighbor
      }
      neighbor.lastVisited = -n
      if(o === 0) {
        toVisit.push(neighbor)
      }
    }
  }
  return null
}

proto.walk = function(point, random) {
  //Alias local properties
  var n = this.vertices.length - 1
  var d = this.dimension
  var verts = this.vertices
  var tuple = this.tuple

  //Compute initial jump cell
  var initIndex = random ? (this.interior.length * Math.random())|0 : (this.interior.length-1)
  var cell = this.interior[ initIndex ]

  //Start walking
outerLoop:
  while(!cell.boundary) {
    var cellVerts = cell.vertices
    var cellAdj = cell.adjacent

    for(var i=0; i<=d; ++i) {
      tuple[i] = verts[cellVerts[i]]
    }
    cell.lastVisited = n

    //Find farthest adjacent cell
    for(var i=0; i<=d; ++i) {
      var neighbor = cellAdj[i]
      if(neighbor.lastVisited >= n) {
        continue
      }
      var prev = tuple[i]
      tuple[i] = point
      var o = this.orient()
      tuple[i] = prev
      if(o < 0) {
        cell = neighbor
        continue outerLoop
      } else {
        if(!neighbor.boundary) {
          neighbor.lastVisited = n
        } else {
          neighbor.lastVisited = -n
        }
      }
    }
    return
  }

  return cell
}

proto.addPeaks = function(point, cell) {
  var n = this.vertices.length - 1
  var d = this.dimension
  var verts = this.vertices
  var tuple = this.tuple
  var interior = this.interior
  var simplices = this.simplices

  //Walking finished at boundary, time to add peaks
  var tovisit = [ cell ]

  //Stretch initial boundary cell into a peak
  cell.lastVisited = n
  cell.vertices[cell.vertices.indexOf(-1)] = n
  cell.boundary = false
  interior.push(cell)

  //Record a list of all new boundaries created by added peaks so we can glue them together when we are all done
  var glueFacets = []

  //Do a traversal of the boundary walking outward from starting peak
  while(tovisit.length > 0) {
    //Pop off peak and walk over adjacent cells
    var cell = tovisit.pop()
    var cellVerts = cell.vertices
    var cellAdj = cell.adjacent
    var indexOfN = cellVerts.indexOf(n)
    if(indexOfN < 0) {
      continue
    }

    for(var i=0; i<=d; ++i) {
      if(i === indexOfN) {
        continue
      }

      //For each boundary neighbor of the cell
      var neighbor = cellAdj[i]
      if(!neighbor.boundary || neighbor.lastVisited >= n) {
        continue
      }

      var nv = neighbor.vertices

      //Test if neighbor is a peak
      if(neighbor.lastVisited !== -n) {      
        //Compute orientation of p relative to each boundary peak
        var indexOfNeg1 = 0
        for(var j=0; j<=d; ++j) {
          if(nv[j] < 0) {
            indexOfNeg1 = j
            tuple[j] = point
          } else {
            tuple[j] = verts[nv[j]]
          }
        }
        var o = this.orient()

        //Test if neighbor cell is also a peak
        if(o > 0) {
          nv[indexOfNeg1] = n
          neighbor.boundary = false
          interior.push(neighbor)
          tovisit.push(neighbor)
          neighbor.lastVisited = n
          continue
        } else {
          neighbor.lastVisited = -n
        }
      }

      var na = neighbor.adjacent

      //Otherwise, replace neighbor with new face
      var vverts = cellVerts.slice()
      var vadj = cellAdj.slice()
      var ncell = new Simplex(vverts, vadj, true)
      simplices.push(ncell)

      //Connect to neighbor
      var opposite = na.indexOf(cell)
      if(opposite < 0) {
        continue
      }
      na[opposite] = ncell
      vadj[indexOfN] = neighbor

      //Connect to cell
      vverts[i] = -1
      vadj[i] = cell
      cellAdj[i] = ncell

      //Flip facet
      ncell.flip()

      //Add to glue list
      for(var j=0; j<=d; ++j) {
        var uu = vverts[j]
        if(uu < 0 || uu === n) {
          continue
        }
        var nface = new Array(d-1)
        var nptr = 0
        for(var k=0; k<=d; ++k) {
          var vv = vverts[k]
          if(vv < 0 || k === j) {
            continue
          }
          nface[nptr++] = vv
        }
        glueFacets.push(new GlueFacet(nface, ncell, j))
      }
    }
  }

  //Glue boundary facets together
  glueFacets.sort(compareGlue)

  for(var i=0; i+1<glueFacets.length; i+=2) {
    var a = glueFacets[i]
    var b = glueFacets[i+1]
    var ai = a.index
    var bi = b.index
    if(ai < 0 || bi < 0) {
      continue
    }
    a.cell.adjacent[a.index] = b.cell
    b.cell.adjacent[b.index] = a.cell
  }
}

proto.insert = function(point, random) {
  //Add point
  var verts = this.vertices
  verts.push(point)

  var cell = this.walk(point, random)
  if(!cell) {
    return
  }

  //Alias local properties
  var d = this.dimension
  var tuple = this.tuple

  //Degenerate case: If point is coplanar to cell, then walk until we find a non-degenerate boundary
  for(var i=0; i<=d; ++i) {
    var vv = cell.vertices[i]
    if(vv < 0) {
      tuple[i] = point
    } else {
      tuple[i] = verts[vv]
    }
  }
  var o = this.orient(tuple)
  if(o < 0) {
    return
  } else if(o === 0) {
    cell = this.handleBoundaryDegeneracy(cell, point)
    if(!cell) {
      return
    }
  }

  //Add peaks
  this.addPeaks(point, cell)
}

//Extract all boundary cells
proto.boundary = function() {
  var d = this.dimension
  var boundary = []
  var cells = this.simplices
  var nc = cells.length
  for(var i=0; i<nc; ++i) {
    var c = cells[i]
    if(c.boundary) {
      var bcell = new Array(d)
      var cv = c.vertices
      var ptr = 0
      var parity = 0
      for(var j=0; j<=d; ++j) {
        if(cv[j] >= 0) {
          bcell[ptr++] = cv[j]
        } else {
          parity = j&1
        }
      }
      if(parity === (d&1)) {
        var t = bcell[0]
        bcell[0] = bcell[1]
        bcell[1] = t
      }
      boundary.push(bcell)
    }
  }
  return boundary
}

function incrementalConvexHull(points, randomSearch) {
  var n = points.length
  if(n === 0) {
    throw new Error("Must have at least d+1 points")
  }
  var d = points[0].length
  if(n <= d) {
    throw new Error("Must input at least d+1 points")
  }

  //FIXME: This could be degenerate, but need to select d+1 non-coplanar points to bootstrap process
  var initialSimplex = points.slice(0, d+1)

  //Make sure initial simplex is positively oriented
  var o = orient.apply(void 0, initialSimplex)
  if(o === 0) {
    throw new Error("Input not in general position")
  }
  var initialCoords = new Array(d+1)
  for(var i=0; i<=d; ++i) {
    initialCoords[i] = i
  }
  if(o < 0) {
    initialCoords[0] = 1
    initialCoords[1] = 0
  }

  //Create initial topological index, glue pointers together (kind of messy)
  var initialCell = new Simplex(initialCoords, new Array(d+1), false)
  var boundary = initialCell.adjacent
  var list = new Array(d+2)
  for(var i=0; i<=d; ++i) {
    var verts = initialCoords.slice()
    for(var j=0; j<=d; ++j) {
      if(j === i) {
        verts[j] = -1
      }
    }
    var t = verts[0]
    verts[0] = verts[1]
    verts[1] = t
    var cell = new Simplex(verts, new Array(d+1), true)
    boundary[i] = cell
    list[i] = cell
  }
  list[d+1] = initialCell
  for(var i=0; i<=d; ++i) {
    var verts = boundary[i].vertices
    var adj = boundary[i].adjacent
    for(var j=0; j<=d; ++j) {
      var v = verts[j]
      if(v < 0) {
        adj[j] = initialCell
        continue
      }
      for(var k=0; k<=d; ++k) {
        if(boundary[k].vertices.indexOf(v) < 0) {
          adj[j] = boundary[k]
        }
      }
    }
  }

  //Initialize triangles
  var triangles = new Triangulation(d, initialSimplex, list)

  //Insert remaining points
  var useRandom = !!randomSearch
  for(var i=d+1; i<n; ++i) {
    triangles.insert(points[i], useRandom)
  }
  
  //Extract boundary cells
  return triangles.boundary()
}
},{"robust-orientation":"/Users/vallette/ANTS/AntColony/node_modules/delaunay-triangulate/node_modules/incremental-convex-hull/node_modules/robust-orientation/orientation.js","simplicial-complex":"/Users/vallette/ANTS/AntColony/node_modules/delaunay-triangulate/node_modules/incremental-convex-hull/node_modules/simplicial-complex/topology.js"}],"/Users/vallette/ANTS/AntColony/node_modules/delaunay-triangulate/node_modules/incremental-convex-hull/node_modules/robust-orientation/node_modules/robust-scale/node_modules/two-sum/two-sum.js":[function(require,module,exports){
"use strict"

module.exports = fastTwoSum

function fastTwoSum(a, b, result) {
	var x = a + b
	var bv = x - a
	var av = x - bv
	var br = b - bv
	var ar = a - av
	if(result) {
		result[0] = ar + br
		result[1] = x
		return result
	}
	return [ar+br, x]
}
},{}],"/Users/vallette/ANTS/AntColony/node_modules/delaunay-triangulate/node_modules/incremental-convex-hull/node_modules/robust-orientation/node_modules/robust-scale/robust-scale.js":[function(require,module,exports){
"use strict"

var twoProduct = require("two-product")
var twoSum = require("two-sum")

module.exports = scaleLinearExpansion

function scaleLinearExpansion(e, scale) {
  var n = e.length
  if(n === 1) {
    var ts = twoProduct(e[0], scale)
    if(ts[0]) {
      return ts
    }
    return [ ts[1] ]
  }
  var g = new Array(2 * n)
  var q = [0.1, 0.1]
  var t = [0.1, 0.1]
  var count = 0
  twoProduct(e[0], scale, q)
  if(q[0]) {
    g[count++] = q[0]
  }
  for(var i=1; i<n; ++i) {
    twoProduct(e[i], scale, t)
    var pq = q[1]
    twoSum(pq, t[0], q)
    if(q[0]) {
      g[count++] = q[0]
    }
    var a = t[1]
    var b = q[1]
    var x = a + b
    var bv = x - a
    var y = b - bv
    q[1] = x
    if(y) {
      g[count++] = y
    }
  }
  if(q[1]) {
    g[count++] = q[1]
  }
  if(count === 0) {
    g[count++] = 0.0
  }
  g.length = count
  return g
}
},{"two-product":"/Users/vallette/ANTS/AntColony/node_modules/delaunay-triangulate/node_modules/incremental-convex-hull/node_modules/robust-orientation/node_modules/two-product/two-product.js","two-sum":"/Users/vallette/ANTS/AntColony/node_modules/delaunay-triangulate/node_modules/incremental-convex-hull/node_modules/robust-orientation/node_modules/robust-scale/node_modules/two-sum/two-sum.js"}],"/Users/vallette/ANTS/AntColony/node_modules/delaunay-triangulate/node_modules/incremental-convex-hull/node_modules/robust-orientation/node_modules/robust-subtract/robust-diff.js":[function(require,module,exports){
"use strict"

module.exports = robustSubtract

//Easy case: Add two scalars
function scalarScalar(a, b) {
  var x = a + b
  var bv = x - a
  var av = x - bv
  var br = b - bv
  var ar = a - av
  var y = ar + br
  if(y) {
    return [y, x]
  }
  return [x]
}

function robustSubtract(e, f) {
  var ne = e.length|0
  var nf = f.length|0
  if(ne === 1 && nf === 1) {
    return scalarScalar(e[0], -f[0])
  }
  var n = ne + nf
  var g = new Array(n)
  var count = 0
  var eptr = 0
  var fptr = 0
  var abs = Math.abs
  var ei = e[eptr]
  var ea = abs(ei)
  var fi = -f[fptr]
  var fa = abs(fi)
  var a, b
  if(ea < fa) {
    b = ei
    eptr += 1
    if(eptr < ne) {
      ei = e[eptr]
      ea = abs(ei)
    }
  } else {
    b = fi
    fptr += 1
    if(fptr < nf) {
      fi = -f[fptr]
      fa = abs(fi)
    }
  }
  if((eptr < ne && ea < fa) || (fptr >= nf)) {
    a = ei
    eptr += 1
    if(eptr < ne) {
      ei = e[eptr]
      ea = abs(ei)
    }
  } else {
    a = fi
    fptr += 1
    if(fptr < nf) {
      fi = -f[fptr]
      fa = abs(fi)
    }
  }
  var x = a + b
  var bv = x - a
  var y = b - bv
  var q0 = y
  var q1 = x
  var _x, _bv, _av, _br, _ar
  while(eptr < ne && fptr < nf) {
    if(ea < fa) {
      a = ei
      eptr += 1
      if(eptr < ne) {
        ei = e[eptr]
        ea = abs(ei)
      }
    } else {
      a = fi
      fptr += 1
      if(fptr < nf) {
        fi = -f[fptr]
        fa = abs(fi)
      }
    }
    b = q0
    x = a + b
    bv = x - a
    y = b - bv
    if(y) {
      g[count++] = y
    }
    _x = q1 + x
    _bv = _x - q1
    _av = _x - _bv
    _br = x - _bv
    _ar = q1 - _av
    q0 = _ar + _br
    q1 = _x
  }
  while(eptr < ne) {
    a = ei
    b = q0
    x = a + b
    bv = x - a
    y = b - bv
    if(y) {
      g[count++] = y
    }
    _x = q1 + x
    _bv = _x - q1
    _av = _x - _bv
    _br = x - _bv
    _ar = q1 - _av
    q0 = _ar + _br
    q1 = _x
    eptr += 1
    if(eptr < ne) {
      ei = e[eptr]
    }
  }
  while(fptr < nf) {
    a = fi
    b = q0
    x = a + b
    bv = x - a
    y = b - bv
    if(y) {
      g[count++] = y
    } 
    _x = q1 + x
    _bv = _x - q1
    _av = _x - _bv
    _br = x - _bv
    _ar = q1 - _av
    q0 = _ar + _br
    q1 = _x
    fptr += 1
    if(fptr < nf) {
      fi = -f[fptr]
    }
  }
  if(q0) {
    g[count++] = q0
  }
  if(q1) {
    g[count++] = q1
  }
  if(!count) {
    g[count++] = 0.0  
  }
  g.length = count
  return g
}
},{}],"/Users/vallette/ANTS/AntColony/node_modules/delaunay-triangulate/node_modules/incremental-convex-hull/node_modules/robust-orientation/node_modules/robust-sum/robust-sum.js":[function(require,module,exports){
"use strict"

module.exports = linearExpansionSum

//Easy case: Add two scalars
function scalarScalar(a, b) {
  var x = a + b
  var bv = x - a
  var av = x - bv
  var br = b - bv
  var ar = a - av
  var y = ar + br
  if(y) {
    return [y, x]
  }
  return [x]
}

function linearExpansionSum(e, f) {
  var ne = e.length|0
  var nf = f.length|0
  if(ne === 1 && nf === 1) {
    return scalarScalar(e[0], f[0])
  }
  var n = ne + nf
  var g = new Array(n)
  var count = 0
  var eptr = 0
  var fptr = 0
  var abs = Math.abs
  var ei = e[eptr]
  var ea = abs(ei)
  var fi = f[fptr]
  var fa = abs(fi)
  var a, b
  if(ea < fa) {
    b = ei
    eptr += 1
    if(eptr < ne) {
      ei = e[eptr]
      ea = abs(ei)
    }
  } else {
    b = fi
    fptr += 1
    if(fptr < nf) {
      fi = f[fptr]
      fa = abs(fi)
    }
  }
  if((eptr < ne && ea < fa) || (fptr >= nf)) {
    a = ei
    eptr += 1
    if(eptr < ne) {
      ei = e[eptr]
      ea = abs(ei)
    }
  } else {
    a = fi
    fptr += 1
    if(fptr < nf) {
      fi = f[fptr]
      fa = abs(fi)
    }
  }
  var x = a + b
  var bv = x - a
  var y = b - bv
  var q0 = y
  var q1 = x
  var _x, _bv, _av, _br, _ar
  while(eptr < ne && fptr < nf) {
    if(ea < fa) {
      a = ei
      eptr += 1
      if(eptr < ne) {
        ei = e[eptr]
        ea = abs(ei)
      }
    } else {
      a = fi
      fptr += 1
      if(fptr < nf) {
        fi = f[fptr]
        fa = abs(fi)
      }
    }
    b = q0
    x = a + b
    bv = x - a
    y = b - bv
    if(y) {
      g[count++] = y
    }
    _x = q1 + x
    _bv = _x - q1
    _av = _x - _bv
    _br = x - _bv
    _ar = q1 - _av
    q0 = _ar + _br
    q1 = _x
  }
  while(eptr < ne) {
    a = ei
    b = q0
    x = a + b
    bv = x - a
    y = b - bv
    if(y) {
      g[count++] = y
    }
    _x = q1 + x
    _bv = _x - q1
    _av = _x - _bv
    _br = x - _bv
    _ar = q1 - _av
    q0 = _ar + _br
    q1 = _x
    eptr += 1
    if(eptr < ne) {
      ei = e[eptr]
    }
  }
  while(fptr < nf) {
    a = fi
    b = q0
    x = a + b
    bv = x - a
    y = b - bv
    if(y) {
      g[count++] = y
    } 
    _x = q1 + x
    _bv = _x - q1
    _av = _x - _bv
    _br = x - _bv
    _ar = q1 - _av
    q0 = _ar + _br
    q1 = _x
    fptr += 1
    if(fptr < nf) {
      fi = f[fptr]
    }
  }
  if(q0) {
    g[count++] = q0
  }
  if(q1) {
    g[count++] = q1
  }
  if(!count) {
    g[count++] = 0.0  
  }
  g.length = count
  return g
}
},{}],"/Users/vallette/ANTS/AntColony/node_modules/delaunay-triangulate/node_modules/incremental-convex-hull/node_modules/robust-orientation/node_modules/two-product/two-product.js":[function(require,module,exports){
"use strict"

module.exports = twoProduct

var SPLITTER = +(Math.pow(2, 27) + 1.0)

function twoProduct(a, b, result) {
  var x = a * b

  var c = SPLITTER * a
  var abig = c - a
  var ahi = c - abig
  var alo = a - ahi

  var d = SPLITTER * b
  var bbig = d - b
  var bhi = d - bbig
  var blo = b - bhi

  var err1 = x - (ahi * bhi)
  var err2 = err1 - (alo * bhi)
  var err3 = err2 - (ahi * blo)

  var y = alo * blo - err3

  if(result) {
    result[0] = y
    result[1] = x
    return result
  }

  return [ y, x ]
}
},{}],"/Users/vallette/ANTS/AntColony/node_modules/delaunay-triangulate/node_modules/incremental-convex-hull/node_modules/robust-orientation/orientation.js":[function(require,module,exports){
"use strict"

var twoProduct = require("two-product")
var robustSum = require("robust-sum")
var robustScale = require("robust-scale")
var robustSubtract = require("robust-subtract")

var NUM_EXPAND = 5

var EPSILON     = 1.1102230246251565e-16
var ERRBOUND3   = (3.0 + 16.0 * EPSILON) * EPSILON
var ERRBOUND4   = (7.0 + 56.0 * EPSILON) * EPSILON

function cofactor(m, c) {
  var result = new Array(m.length-1)
  for(var i=1; i<m.length; ++i) {
    var r = result[i-1] = new Array(m.length-1)
    for(var j=0,k=0; j<m.length; ++j) {
      if(j === c) {
        continue
      }
      r[k++] = m[i][j]
    }
  }
  return result
}

function matrix(n) {
  var result = new Array(n)
  for(var i=0; i<n; ++i) {
    result[i] = new Array(n)
    for(var j=0; j<n; ++j) {
      result[i][j] = ["m", j, "[", (n-i-1), "]"].join("")
    }
  }
  return result
}

function sign(n) {
  if(n & 1) {
    return "-"
  }
  return ""
}

function generateSum(expr) {
  if(expr.length === 1) {
    return expr[0]
  } else if(expr.length === 2) {
    return ["sum(", expr[0], ",", expr[1], ")"].join("")
  } else {
    var m = expr.length>>1
    return ["sum(", generateSum(expr.slice(0, m)), ",", generateSum(expr.slice(m)), ")"].join("")
  }
}

function determinant(m) {
  if(m.length === 2) {
    return [["sum(prod(", m[0][0], ",", m[1][1], "),prod(-", m[0][1], ",", m[1][0], "))"].join("")]
  } else {
    var expr = []
    for(var i=0; i<m.length; ++i) {
      expr.push(["scale(", generateSum(determinant(cofactor(m, i))), ",", sign(i), m[0][i], ")"].join(""))
    }
    return expr
  }
}

function orientation(n) {
  var pos = []
  var neg = []
  var m = matrix(n)
  var args = []
  for(var i=0; i<n; ++i) {
    if((i&1)===0) {
      pos.push.apply(pos, determinant(cofactor(m, i)))
    } else {
      neg.push.apply(neg, determinant(cofactor(m, i)))
    }
    args.push("m" + i)
  }
  var posExpr = generateSum(pos)
  var negExpr = generateSum(neg)
  var funcName = "orientation" + n + "Exact"
  var code = ["function ", funcName, "(", args.join(), "){var p=", posExpr, ",n=", negExpr, ",d=sub(p,n);\
return d[d.length-1];};return ", funcName].join("")
  var proc = new Function("sum", "prod", "scale", "sub", code)
  return proc(robustSum, twoProduct, robustScale, robustSubtract)
}

var orientation3Exact = orientation(3)
var orientation4Exact = orientation(4)

var CACHED = [
  function orientation0() { return 0 },
  function orientation1() { return 0 },
  function orientation2(a, b) { 
    return b[0] - a[0]
  },
  function orientation3(a, b, c) {
    var l = (a[1] - c[1]) * (b[0] - c[0])
    var r = (a[0] - c[0]) * (b[1] - c[1])
    var det = l - r
    var s
    if(l > 0) {
      if(r <= 0) {
        return det
      } else {
        s = l + r
      }
    } else if(l < 0) {
      if(r >= 0) {
        return det
      } else {
        s = -(l + r)
      }
    } else {
      return det
    }
    var tol = ERRBOUND3 * s
    if(det >= tol || det <= -tol) {
      return det
    }
    return orientation3Exact(a, b, c)
  },
  function orientation4(a,b,c,d) {
    var adx = a[0] - d[0]
    var bdx = b[0] - d[0]
    var cdx = c[0] - d[0]
    var ady = a[1] - d[1]
    var bdy = b[1] - d[1]
    var cdy = c[1] - d[1]
    var adz = a[2] - d[2]
    var bdz = b[2] - d[2]
    var cdz = c[2] - d[2]
    var bdxcdy = bdx * cdy
    var cdxbdy = cdx * bdy
    var cdxady = cdx * ady
    var adxcdy = adx * cdy
    var adxbdy = adx * bdy
    var bdxady = bdx * ady
    var det = adz * (bdxcdy - cdxbdy) 
            + bdz * (cdxady - adxcdy)
            + cdz * (adxbdy - bdxady)
    var permanent = (Math.abs(bdxcdy) + Math.abs(cdxbdy)) * Math.abs(adz)
                  + (Math.abs(cdxady) + Math.abs(adxcdy)) * Math.abs(bdz)
                  + (Math.abs(adxbdy) + Math.abs(bdxady)) * Math.abs(cdz)
    var tol = ERRBOUND4 * permanent
    if ((det > tol) || (-det > tol)) {
      return det
    }
    return orientation4Exact(a,b,c,d)
  }
]

function slowOrient(args) {
  var proc = CACHED[args.length]
  if(!proc) {
    proc = CACHED[args.length] = orientation(args.length)
  }
  return proc.apply(undefined, args)
}

function generateOrientationProc() {
  while(CACHED.length <= NUM_EXPAND) {
    CACHED.push(orientation(CACHED.length))
  }
  var args = []
  var procArgs = ["slow"]
  for(var i=0; i<=NUM_EXPAND; ++i) {
    args.push("a" + i)
    procArgs.push("o" + i)
  }
  var code = [
    "function getOrientation(", args.join(), "){switch(arguments.length){case 0:case 1:return 0;"
  ]
  for(var i=2; i<=NUM_EXPAND; ++i) {
    code.push("case ", i, ":return o", i, "(", args.slice(0, i).join(), ");")
  }
  code.push("}var s=new Array(arguments.length);for(var i=0;i<arguments.length;++i){s[i]=arguments[i]};return slow(s);}return getOrientation")
  procArgs.push(code.join(""))

  var proc = Function.apply(undefined, procArgs)
  module.exports = proc.apply(undefined, [slowOrient].concat(CACHED))
  for(var i=0; i<=NUM_EXPAND; ++i) {
    module.exports[i] = CACHED[i]
  }
}

generateOrientationProc()
},{"robust-scale":"/Users/vallette/ANTS/AntColony/node_modules/delaunay-triangulate/node_modules/incremental-convex-hull/node_modules/robust-orientation/node_modules/robust-scale/robust-scale.js","robust-subtract":"/Users/vallette/ANTS/AntColony/node_modules/delaunay-triangulate/node_modules/incremental-convex-hull/node_modules/robust-orientation/node_modules/robust-subtract/robust-diff.js","robust-sum":"/Users/vallette/ANTS/AntColony/node_modules/delaunay-triangulate/node_modules/incremental-convex-hull/node_modules/robust-orientation/node_modules/robust-sum/robust-sum.js","two-product":"/Users/vallette/ANTS/AntColony/node_modules/delaunay-triangulate/node_modules/incremental-convex-hull/node_modules/robust-orientation/node_modules/two-product/two-product.js"}],"/Users/vallette/ANTS/AntColony/node_modules/delaunay-triangulate/node_modules/incremental-convex-hull/node_modules/simplicial-complex/node_modules/bit-twiddle/twiddle.js":[function(require,module,exports){
/**
 * Bit twiddling hacks for JavaScript.
 *
 * Author: Mikola Lysenko
 *
 * Ported from Stanford bit twiddling hack library:
 *    http://graphics.stanford.edu/~seander/bithacks.html
 */

"use strict"; "use restrict";

//Number of bits in an integer
var INT_BITS = 32;

//Constants
exports.INT_BITS  = INT_BITS;
exports.INT_MAX   =  0x7fffffff;
exports.INT_MIN   = -1<<(INT_BITS-1);

//Returns -1, 0, +1 depending on sign of x
exports.sign = function(v) {
  return (v > 0) - (v < 0);
}

//Computes absolute value of integer
exports.abs = function(v) {
  var mask = v >> (INT_BITS-1);
  return (v ^ mask) - mask;
}

//Computes minimum of integers x and y
exports.min = function(x, y) {
  return y ^ ((x ^ y) & -(x < y));
}

//Computes maximum of integers x and y
exports.max = function(x, y) {
  return x ^ ((x ^ y) & -(x < y));
}

//Checks if a number is a power of two
exports.isPow2 = function(v) {
  return !(v & (v-1)) && (!!v);
}

//Computes log base 2 of v
exports.log2 = function(v) {
  var r, shift;
  r =     (v > 0xFFFF) << 4; v >>>= r;
  shift = (v > 0xFF  ) << 3; v >>>= shift; r |= shift;
  shift = (v > 0xF   ) << 2; v >>>= shift; r |= shift;
  shift = (v > 0x3   ) << 1; v >>>= shift; r |= shift;
  return r | (v >> 1);
}

//Computes log base 10 of v
exports.log10 = function(v) {
  return  (v >= 1000000000) ? 9 : (v >= 100000000) ? 8 : (v >= 10000000) ? 7 :
          (v >= 1000000) ? 6 : (v >= 100000) ? 5 : (v >= 10000) ? 4 :
          (v >= 1000) ? 3 : (v >= 100) ? 2 : (v >= 10) ? 1 : 0;
}

//Counts number of bits
exports.popCount = function(v) {
  v = v - ((v >>> 1) & 0x55555555);
  v = (v & 0x33333333) + ((v >>> 2) & 0x33333333);
  return ((v + (v >>> 4) & 0xF0F0F0F) * 0x1010101) >>> 24;
}

//Counts number of trailing zeros
function countTrailingZeros(v) {
  var c = 32;
  v &= -v;
  if (v) c--;
  if (v & 0x0000FFFF) c -= 16;
  if (v & 0x00FF00FF) c -= 8;
  if (v & 0x0F0F0F0F) c -= 4;
  if (v & 0x33333333) c -= 2;
  if (v & 0x55555555) c -= 1;
  return c;
}
exports.countTrailingZeros = countTrailingZeros;

//Rounds to next power of 2
exports.nextPow2 = function(v) {
  v += v === 0;
  --v;
  v |= v >>> 1;
  v |= v >>> 2;
  v |= v >>> 4;
  v |= v >>> 8;
  v |= v >>> 16;
  return v + 1;
}

//Rounds down to previous power of 2
exports.prevPow2 = function(v) {
  v |= v >>> 1;
  v |= v >>> 2;
  v |= v >>> 4;
  v |= v >>> 8;
  v |= v >>> 16;
  return v - (v>>>1);
}

//Computes parity of word
exports.parity = function(v) {
  v ^= v >>> 16;
  v ^= v >>> 8;
  v ^= v >>> 4;
  v &= 0xf;
  return (0x6996 >>> v) & 1;
}

var REVERSE_TABLE = new Array(256);

(function(tab) {
  for(var i=0; i<256; ++i) {
    var v = i, r = i, s = 7;
    for (v >>>= 1; v; v >>>= 1) {
      r <<= 1;
      r |= v & 1;
      --s;
    }
    tab[i] = (r << s) & 0xff;
  }
})(REVERSE_TABLE);

//Reverse bits in a 32 bit word
exports.reverse = function(v) {
  return  (REVERSE_TABLE[ v         & 0xff] << 24) |
          (REVERSE_TABLE[(v >>> 8)  & 0xff] << 16) |
          (REVERSE_TABLE[(v >>> 16) & 0xff] << 8)  |
           REVERSE_TABLE[(v >>> 24) & 0xff];
}

//Interleave bits of 2 coordinates with 16 bits.  Useful for fast quadtree codes
exports.interleave2 = function(x, y) {
  x &= 0xFFFF;
  x = (x | (x << 8)) & 0x00FF00FF;
  x = (x | (x << 4)) & 0x0F0F0F0F;
  x = (x | (x << 2)) & 0x33333333;
  x = (x | (x << 1)) & 0x55555555;

  y &= 0xFFFF;
  y = (y | (y << 8)) & 0x00FF00FF;
  y = (y | (y << 4)) & 0x0F0F0F0F;
  y = (y | (y << 2)) & 0x33333333;
  y = (y | (y << 1)) & 0x55555555;

  return x | (y << 1);
}

//Extracts the nth interleaved component
exports.deinterleave2 = function(v, n) {
  v = (v >>> n) & 0x55555555;
  v = (v | (v >>> 1))  & 0x33333333;
  v = (v | (v >>> 2))  & 0x0F0F0F0F;
  v = (v | (v >>> 4))  & 0x00FF00FF;
  v = (v | (v >>> 16)) & 0x000FFFF;
  return (v << 16) >> 16;
}


//Interleave bits of 3 coordinates, each with 10 bits.  Useful for fast octree codes
exports.interleave3 = function(x, y, z) {
  x &= 0x3FF;
  x  = (x | (x<<16)) & 4278190335;
  x  = (x | (x<<8))  & 251719695;
  x  = (x | (x<<4))  & 3272356035;
  x  = (x | (x<<2))  & 1227133513;

  y &= 0x3FF;
  y  = (y | (y<<16)) & 4278190335;
  y  = (y | (y<<8))  & 251719695;
  y  = (y | (y<<4))  & 3272356035;
  y  = (y | (y<<2))  & 1227133513;
  x |= (y << 1);
  
  z &= 0x3FF;
  z  = (z | (z<<16)) & 4278190335;
  z  = (z | (z<<8))  & 251719695;
  z  = (z | (z<<4))  & 3272356035;
  z  = (z | (z<<2))  & 1227133513;
  
  return x | (z << 2);
}

//Extracts nth interleaved component of a 3-tuple
exports.deinterleave3 = function(v, n) {
  v = (v >>> n)       & 1227133513;
  v = (v | (v>>>2))   & 3272356035;
  v = (v | (v>>>4))   & 251719695;
  v = (v | (v>>>8))   & 4278190335;
  v = (v | (v>>>16))  & 0x3FF;
  return (v<<22)>>22;
}

//Computes next combination in colexicographic order (this is mistakenly called nextPermutation on the bit twiddling hacks page)
exports.nextCombination = function(v) {
  var t = v | (v - 1);
  return (t + 1) | (((~t & -~t) - 1) >>> (countTrailingZeros(v) + 1));
}


},{}],"/Users/vallette/ANTS/AntColony/node_modules/delaunay-triangulate/node_modules/incremental-convex-hull/node_modules/simplicial-complex/node_modules/union-find/index.js":[function(require,module,exports){
"use strict"; "use restrict";

module.exports = UnionFind;

function UnionFind(count) {
  this.roots = new Array(count);
  this.ranks = new Array(count);
  
  for(var i=0; i<count; ++i) {
    this.roots[i] = i;
    this.ranks[i] = 0;
  }
}

var proto = UnionFind.prototype

Object.defineProperty(proto, "length", {
  "get": function() {
    return this.roots.length
  }
})

proto.makeSet = function() {
  var n = this.roots.length;
  this.roots.push(n);
  this.ranks.push(0);
  return n;
}

proto.find = function(x) {
  var roots = this.roots;
  while(roots[x] !== x) {
    var y = roots[x];
    roots[x] = roots[y];
    x = y;
  }
  return x;
}

proto.link = function(x, y) {
  var xr = this.find(x)
    , yr = this.find(y);
  if(xr === yr) {
    return;
  }
  var ranks = this.ranks
    , roots = this.roots
    , xd    = ranks[xr]
    , yd    = ranks[yr];
  if(xd < yd) {
    roots[xr] = yr;
  } else if(yd < xd) {
    roots[yr] = xr;
  } else {
    roots[yr] = xr;
    ++ranks[xr];
  }
}
},{}],"/Users/vallette/ANTS/AntColony/node_modules/delaunay-triangulate/node_modules/incremental-convex-hull/node_modules/simplicial-complex/topology.js":[function(require,module,exports){
"use strict"; "use restrict";

var bits      = require("bit-twiddle")
  , UnionFind = require("union-find")

//Returns the dimension of a cell complex
function dimension(cells) {
  var d = 0
    , max = Math.max
  for(var i=0, il=cells.length; i<il; ++i) {
    d = max(d, cells[i].length)
  }
  return d-1
}
exports.dimension = dimension

//Counts the number of vertices in faces
function countVertices(cells) {
  var vc = -1
    , max = Math.max
  for(var i=0, il=cells.length; i<il; ++i) {
    var c = cells[i]
    for(var j=0, jl=c.length; j<jl; ++j) {
      vc = max(vc, c[j])
    }
  }
  return vc+1
}
exports.countVertices = countVertices

//Returns a deep copy of cells
function cloneCells(cells) {
  var ncells = new Array(cells.length)
  for(var i=0, il=cells.length; i<il; ++i) {
    ncells[i] = cells[i].slice(0)
  }
  return ncells
}
exports.cloneCells = cloneCells

//Ranks a pair of cells up to permutation
function compareCells(a, b) {
  var n = a.length
    , t = a.length - b.length
    , min = Math.min
  if(t) {
    return t
  }
  switch(n) {
    case 0:
      return 0;
    case 1:
      return a[0] - b[0];
    case 2:
      var d = a[0]+a[1]-b[0]-b[1]
      if(d) {
        return d
      }
      return min(a[0],a[1]) - min(b[0],b[1])
    case 3:
      var l1 = a[0]+a[1]
        , m1 = b[0]+b[1]
      d = l1+a[2] - (m1+b[2])
      if(d) {
        return d
      }
      var l0 = min(a[0], a[1])
        , m0 = min(b[0], b[1])
        , d  = min(l0, a[2]) - min(m0, b[2])
      if(d) {
        return d
      }
      return min(l0+a[2], l1) - min(m0+b[2], m1)
    
    //TODO: Maybe optimize n=4 as well?
    
    default:
      var as = a.slice(0)
      as.sort()
      var bs = b.slice(0)
      bs.sort()
      for(var i=0; i<n; ++i) {
        t = as[i] - bs[i]
        if(t) {
          return t
        }
      }
      return 0
  }
}
exports.compareCells = compareCells

function compareZipped(a, b) {
  return compareCells(a[0], b[0])
}

//Puts a cell complex into normal order for the purposes of findCell queries
function normalize(cells, attr) {
  if(attr) {
    var len = cells.length
    var zipped = new Array(len)
    for(var i=0; i<len; ++i) {
      zipped[i] = [cells[i], attr[i]]
    }
    zipped.sort(compareZipped)
    for(var i=0; i<len; ++i) {
      cells[i] = zipped[i][0]
      attr[i] = zipped[i][1]
    }
    return cells
  } else {
    cells.sort(compareCells)
    return cells
  }
}
exports.normalize = normalize

//Removes all duplicate cells in the complex
function unique(cells) {
  if(cells.length === 0) {
    return []
  }
  var ptr = 1
    , len = cells.length
  for(var i=1; i<len; ++i) {
    var a = cells[i]
    if(compareCells(a, cells[i-1])) {
      if(i === ptr) {
        ptr++
        continue
      }
      cells[ptr++] = a
    }
  }
  cells.length = ptr
  return cells
}
exports.unique = unique;

//Finds a cell in a normalized cell complex
function findCell(cells, c) {
  var lo = 0
    , hi = cells.length-1
    , r  = -1
  while (lo <= hi) {
    var mid = (lo + hi) >> 1
      , s   = compareCells(cells[mid], c)
    if(s <= 0) {
      if(s === 0) {
        r = mid
      }
      lo = mid + 1
    } else if(s > 0) {
      hi = mid - 1
    }
  }
  return r
}
exports.findCell = findCell;

//Builds an index for an n-cell.  This is more general than dual, but less efficient
function incidence(from_cells, to_cells) {
  var index = new Array(from_cells.length)
  for(var i=0, il=index.length; i<il; ++i) {
    index[i] = []
  }
  var b = []
  for(var i=0, n=to_cells.length; i<n; ++i) {
    var c = to_cells[i]
    var cl = c.length
    for(var k=1, kn=(1<<cl); k<kn; ++k) {
      b.length = bits.popCount(k)
      var l = 0
      for(var j=0; j<cl; ++j) {
        if(k & (1<<j)) {
          b[l++] = c[j]
        }
      }
      var idx=findCell(from_cells, b)
      if(idx < 0) {
        continue
      }
      while(true) {
        index[idx++].push(i)
        if(idx >= from_cells.length || compareCells(from_cells[idx], b) !== 0) {
          break
        }
      }
    }
  }
  return index
}
exports.incidence = incidence

//Computes the dual of the mesh.  This is basically an optimized version of buildIndex for the situation where from_cells is just the list of vertices
function dual(cells, vertex_count) {
  if(!vertex_count) {
    return incidence(unique(skeleton(cells, 0)), cells, 0)
  }
  var res = new Array(vertex_count)
  for(var i=0; i<vertex_count; ++i) {
    res[i] = []
  }
  for(var i=0, len=cells.length; i<len; ++i) {
    var c = cells[i]
    for(var j=0, cl=c.length; j<cl; ++j) {
      res[c[j]].push(i)
    }
  }
  return res
}
exports.dual = dual

//Enumerates all cells in the complex
function explode(cells) {
  var result = []
  for(var i=0, il=cells.length; i<il; ++i) {
    var c = cells[i]
      , cl = c.length|0
    for(var j=1, jl=(1<<cl); j<jl; ++j) {
      var b = []
      for(var k=0; k<cl; ++k) {
        if((j >>> k) & 1) {
          b.push(c[k])
        }
      }
      result.push(b)
    }
  }
  return normalize(result)
}
exports.explode = explode

//Enumerates all of the n-cells of a cell complex
function skeleton(cells, n) {
  if(n < 0) {
    return []
  }
  var result = []
    , k0     = (1<<(n+1))-1
  for(var i=0; i<cells.length; ++i) {
    var c = cells[i]
    for(var k=k0; k<(1<<c.length); k=bits.nextCombination(k)) {
      var b = new Array(n+1)
        , l = 0
      for(var j=0; j<c.length; ++j) {
        if(k & (1<<j)) {
          b[l++] = c[j]
        }
      }
      result.push(b)
    }
  }
  return normalize(result)
}
exports.skeleton = skeleton;

//Computes the boundary of all cells, does not remove duplicates
function boundary(cells) {
  var res = []
  for(var i=0,il=cells.length; i<il; ++i) {
    var c = cells[i]
    for(var j=0,cl=c.length; j<cl; ++j) {
      var b = new Array(c.length-1)
      for(var k=0, l=0; k<cl; ++k) {
        if(k !== j) {
          b[l++] = c[k]
        }
      }
      res.push(b)
    }
  }
  return normalize(res)
}
exports.boundary = boundary;

//Computes connected components for a dense cell complex
function connectedComponents_dense(cells, vertex_count) {
  var labels = new UnionFind(vertex_count)
  for(var i=0; i<cells.length; ++i) {
    var c = cells[i]
    for(var j=0; j<c.length; ++j) {
      for(var k=j+1; k<c.length; ++k) {
        labels.link(c[j], c[k])
      }
    }
  }
  var components = []
    , component_labels = labels.ranks
  for(var i=0; i<component_labels.length; ++i) {
    component_labels[i] = -1
  }
  for(var i=0; i<cells.length; ++i) {
    var l = labels.find(cells[i][0])
    if(component_labels[l] < 0) {
      component_labels[l] = components.length
      components.push([cells[i].slice(0)])
    } else {
      components[component_labels[l]].push(cells[i].slice(0))
    }
  }
  return components
}

//Computes connected components for a sparse graph
function connectedComponents_sparse(cells) {
  var vertices  = unique(normalize(skeleton(cells, 0)))
    , labels    = new UnionFind(vertices.length)
  for(var i=0; i<cells.length; ++i) {
    var c = cells[i]
    for(var j=0; j<c.length; ++j) {
      var vj = findCell(vertices, [c[j]])
      for(var k=j+1; k<c.length; ++k) {
        labels.link(vj, findCell(vertices, [c[k]]))
      }
    }
  }
  var components        = []
    , component_labels  = labels.ranks
  for(var i=0; i<component_labels.length; ++i) {
    component_labels[i] = -1
  }
  for(var i=0; i<cells.length; ++i) {
    var l = labels.find(findCell(vertices, [cells[i][0]]));
    if(component_labels[l] < 0) {
      component_labels[l] = components.length
      components.push([cells[i].slice(0)])
    } else {
      components[component_labels[l]].push(cells[i].slice(0))
    }
  }
  return components
}

//Computes connected components for a cell complex
function connectedComponents(cells, vertex_count) {
  if(vertex_count) {
    return connectedComponents_dense(cells, vertex_count)
  }
  return connectedComponents_sparse(cells)
}
exports.connectedComponents = connectedComponents

},{"bit-twiddle":"/Users/vallette/ANTS/AntColony/node_modules/delaunay-triangulate/node_modules/incremental-convex-hull/node_modules/simplicial-complex/node_modules/bit-twiddle/twiddle.js","union-find":"/Users/vallette/ANTS/AntColony/node_modules/delaunay-triangulate/node_modules/incremental-convex-hull/node_modules/simplicial-complex/node_modules/union-find/index.js"}],"/Users/vallette/ANTS/AntColony/node_modules/delaunay-triangulate/node_modules/uniq/uniq.js":[function(require,module,exports){
"use strict"

function unique_pred(list, compare) {
  var ptr = 1
    , len = list.length
    , a=list[0], b=list[0]
  for(var i=1; i<len; ++i) {
    b = a
    a = list[i]
    if(compare(a, b)) {
      if(i === ptr) {
        ptr++
        continue
      }
      list[ptr++] = a
    }
  }
  list.length = ptr
  return list
}

function unique_eq(list) {
  var ptr = 1
    , len = list.length
    , a=list[0], b = list[0]
  for(var i=1; i<len; ++i, b=a) {
    b = a
    a = list[i]
    if(a !== b) {
      if(i === ptr) {
        ptr++
        continue
      }
      list[ptr++] = a
    }
  }
  list.length = ptr
  return list
}

function unique(list, compare, sorted) {
  if(list.length === 0) {
    return list
  }
  if(compare) {
    if(!sorted) {
      list.sort(compare)
    }
    return unique_pred(list, compare)
  }
  if(!sorted) {
    list.sort()
  }
  return unique_eq(list)
}

module.exports = unique

},{}],"/Users/vallette/ANTS/AntColony/node_modules/delaunay-triangulate/triangulate.js":[function(require,module,exports){
"use strict"

var ch = require("incremental-convex-hull")
var uniq = require("uniq")

module.exports = triangulate

function LiftedPoint(p, i) {
  this.point = p
  this.index = i
}

function compareLifted(a, b) {
  var ap = a.point
  var bp = b.point
  var d = ap.length
  for(var i=0; i<d; ++i) {
    var s = bp[i] - ap[i]
    if(s) {
      return s
    }
  }
  return 0
}

function triangulate1D(n, points, includePointAtInfinity) {
  if(n === 1) {
    if(includePointAtInfinity) {
      return [ [-1, 0] ]
    } else {
      return []
    }
  }
  var lifted = points.map(function(p, i) {
    return [ p[0], i ]
  })
  lifted.sort(function(a,b) {
    return a[0] - b[0]
  })
  var cells = new Array(n - 1)
  for(var i=1; i<n; ++i) {
    var a = lifted[i-1]
    var b = lifted[i]
    cells[i-1] = [ a[1], b[1] ]
  }
  if(includePointAtInfinity) {
    cells.push(
      [ -1, cells[0][1], ],
      [ cells[n-1][1], -1 ])
  }
  return cells
}

function triangulate(points, includePointAtInfinity) {
  var n = points.length
  if(n === 0) {
    return []
  }
  
  var d = points[0].length
  if(d < 1) {
    return []
  }

  //Special case:  For 1D we can just sort the points
  if(d === 1) {
    return triangulate1D(n, points, includePointAtInfinity)
  }
  
  //Lift points, sort
  var lifted = new Array(n)
  var upper = 1.0
  for(var i=0; i<n; ++i) {
    var p = points[i]
    var x = new Array(d+1)
    var l = 0.0
    for(var j=0; j<d; ++j) {
      var v = p[j]
      x[j] = v
      l += v * v
    }
    x[d] = l
    lifted[i] = new LiftedPoint(x, i)
    upper = Math.max(l, upper)
  }
  uniq(lifted, compareLifted)
  
  //Double points
  n = lifted.length

  //Create new list of points
  var dpoints = new Array(n + d + 1)
  var dindex = new Array(n + d + 1)

  //Add steiner points at top
  var u = (d+1) * (d+1) * upper
  var y = new Array(d+1)
  for(var i=0; i<=d; ++i) {
    y[i] = 0.0
  }
  y[d] = u

  dpoints[0] = y.slice()
  dindex[0] = -1

  for(var i=0; i<=d; ++i) {
    var x = y.slice()
    x[i] = 1
    dpoints[i+1] = x
    dindex[i+1] = -1
  }

  //Copy rest of the points over
  for(var i=0; i<n; ++i) {
    var h = lifted[i]
    dpoints[i + d + 1] = h.point
    dindex[i + d + 1] =  h.index
  }

  //Construct convex hull
  var hull = ch(dpoints, false)
  if(includePointAtInfinity) {
    hull = hull.filter(function(cell) {
      var count = 0
      for(var j=0; j<=d; ++j) {
        var v = dindex[cell[j]]
        if(v < 0) {
          if(++count >= 2) {
            return false
          }
        }
        cell[j] = v
      }
      return true
    })
  } else {
    hull = hull.filter(function(cell) {
      for(var i=0; i<=d; ++i) {
        var v = dindex[cell[i]]
        if(v < 0) {
          return false
        }
        cell[i] = v
      }
      return true
    })
  }

  if(d & 1) {
    for(var i=0; i<hull.length; ++i) {
      var h = hull[i]
      var x = h[0]
      h[0] = h[1]
      h[1] = x
    }
  }

  return hull
}
},{"incremental-convex-hull":"/Users/vallette/ANTS/AntColony/node_modules/delaunay-triangulate/node_modules/incremental-convex-hull/ich.js","uniq":"/Users/vallette/ANTS/AntColony/node_modules/delaunay-triangulate/node_modules/uniq/uniq.js"}],"/Users/vallette/ANTS/AntColony/node_modules/parse-svg-path/index.js":[function(require,module,exports){

module.exports = parse

/**
 * expected argument lengths
 * @type {Object}
 */

var length = {a: 7, c: 6, h: 1, l: 2, m: 2, q: 4, s: 4, t: 2, v: 1, z: 0}

/**
 * segment pattern
 * @type {RegExp}
 */

var segment = /([astvzqmhlc])([^astvzqmhlc]*)/ig

/**
 * parse an svg path data string. Generates an Array
 * of commands where each command is an Array of the
 * form `[command, arg1, arg2, ...]`
 *
 * @param {String} path
 * @return {Array}
 */

function parse(path) {
	var data = []
	path.replace(segment, function(_, command, args){
		var type = command.toLowerCase()
		args = parseValues(args)

		// overloaded moveTo
		if (type == 'm' && args.length > 2) {
			data.push([command].concat(args.splice(0, 2)))
			type = 'l'
			command = command == 'm' ? 'l' : 'L'
		}

		while (true) {
			if (args.length == length[type]) {
				args.unshift(command)
				return data.push(args)
			}
			if (args.length < length[type]) throw new Error('malformed path data')
			data.push([command].concat(args.splice(0, length[type])))
		}
	})
	return data
}

function parseValues(args){
	args = args.match(/-?[.0-9]+(?:e[-+]?\d+)?/ig)
	return args ? args.map(Number) : []
}

},{}],"/Users/vallette/ANTS/AntColony/src/ant.js":[function(require,module,exports){
'use strict';

var sign = require('./utilities.js').sign;
var calculateDistance = require('./utilities.js').distance;

var points = require('./initializePoints.js').points;
var citySet = require('./initializePoints.js').citySet;
var textPointsId = require('./initializePoints.js').textPointsId;
var possibleStartPointsId = require('./initializePoints.js').possibleStartPointsId;

var liveMousePosition = require('./mouse.js');

var Vector = require('./vector.js');

var random = Math.random;
var floor = Math.floor;
var REPULSION = 0.05;


module.exports = function(container){

    var mouse = liveMousePosition(container);


    function Ant(point) {
        this.x = point.x;                
        this.y = point.y;
        this.velocity = 0.003;
        this.edge = undefined;
        this.state = "forage";
        this.edges = [];
        this.lastCity = undefined;
        this.origin = point;
        this.destination = undefined;
        this.orientation = undefined;
        this.direction = new Vector(0,0);
        this.prog = 0;
    }
    // forage: the ant wanders around without any pheromon deposition
    // once it finds a city, it starts remembering the nodes it goes through
    // when it finds another city, it computes the path length and adds pheromons one each edges
    // proportionnaly to the shortestness of the path
    // it resets the list of nodes and continues
    // while foraging the ant choses the path with a pheromon preference


    // static methods
    Ant.generateRandStartPoint = function() {
        var randId = Math.floor(possibleStartPointsId.length * random());
        var randStartPoint = points[possibleStartPointsId[randId]];
        return randStartPoint;
    }


    // methods
    Ant.prototype = {

        transit: function(){
            switch (this.state) {
            case "forage":
                var res = this.move();
                if (res.cityReached) {
                    this.state = "pheromoning";
                    this.lastCity = this.origin.id;
                };
                break;
            case "pheromoning":
                var res = this.move();
                if (res.edgeChanged) {
                    this.edges.push(this.edge);
                    // found a city
                    if (res.cityReached && (this.origin.id != this.lastCity) ){
                        // compute the length of the path
                        var pathLength = this.edges.map(function(e){return e.distance}).reduce(function(a,b){return a + b});
                        var deltaPheromone = 1/pathLength;
                        this.edges.forEach(function(e){
                            var a = e.pt1, b = e.pt2, weight = 1;  
                            // increased dropped pheromons for textEdges
                            if ((citySet.indexOf(a.id) != -1) && citySet.indexOf(b.id) != -1 && (Math.abs(a.id - b.id) == 1))
                            {
                                weight *= 10;
                            }
                            e.pheromon += (deltaPheromone * weight);
                        });

                        this.edges = [this.edge];
                        this.lastCity = this.origin.id;
                    }
                }
              break;
            }

        },

        setDirection: function(){
            var possibleEdges = [];

            for (var i = 0; i < this.origin.nexts.length; i++)
            {
                possibleEdges[i] = this.origin.nexts[i];
            } 

            possibleEdges.splice(possibleEdges.indexOf(this.edge),1);

            // flip a coin and either take the smelliest path or a random one
            if (random() > 0.5){
                var smells = possibleEdges.map(function(e){return e.pheromon});
                var index = smells.indexOf(Math.max.apply(Math, smells));
                this.edge = possibleEdges[index];
            } 
            else
                this.edge = possibleEdges[floor(random()*possibleEdges.length)];

            // set the destination point, being edge.pt1 or edge.pt2
            this.destination = (this.origin == this.edge.pt1) ? this.edge.pt2 : this.edge.pt1;

            this.direction.x = this.destination.x - this.origin.x; 
            this.direction.y = this.destination.y - this.origin.y;

            this.direction.normalize();
        },

        move: function(){
            var edgeChanged;
            var cityReached = false;

            this.direction.x = this.destination.x - this.x; 
            this.direction.y = this.destination.y - this.y;
            this.direction.normalize();

            // on edge
            // if ((this.prog < this.edge.distance) || (calculateDistance(this, this.destination) > 0.001){
            if ((calculateDistance(this, this.destination) > 0.005)){

                // a delta movement will be applied if collision with obstacle detected
                var delta = this.avoidObstacle();

                this.x += this.velocity * this.direction.x + delta.x * 0.005;
                this.y += this.velocity * this.direction.y + delta.y * 0.005;

                this.prog = this.calculateProgression();
                // this.prog = calculateDistance(this, this.origin);
                //console.log(this.prog / this.edge.distance);
                
                edgeChanged = false;

            // on vertex
            } else {
                this.step = 0;
                this.prog = 0;
                this.origin = this.destination;
                this.x = this.origin.x;
                this.y = this.origin.y;

                this.setDirection();

                cityReached = (citySet.indexOf(this.origin.id) != -1);
                edgeChanged = true;
            }
            return {cityReached: cityReached, edgeChanged: edgeChanged};
        },

        avoidObstacle: function(){
            var distance = calculateDistance(this, mouse);
        
            if (distance <= REPULSION)
            {
                // if (distanceEdge > 0.01){
                //     this.direction.x = this.destination.x - this.x; 
                //     this.direction.y = this.destination.y - this.y;
                //     this.direction.normalize();
                // }

                return {
                    // delta movement is composed of a repulsion delta and a circular delta 
                    x: (this.x - mouse.x)/distance + (this.y - mouse.y)/distance * 1,
                    y: (this.y - mouse.y)/distance - (this.x - mouse.x)/distance * 1
                    // x: (this.x - mouse.x)/distance + (this.y - mouse.y)/distance * 1 + 0.1 * this.direction.x,
                    // y: (this.y - mouse.y)/distance - (this.x - mouse.x)/distance * 1 + 0.1 * this.direction.y 
                };
            }
            else
                return {x:0, y:0};
        },

        calculateProgression: function(){
            var v = new Vector(this.x - this.origin.x, this.y - this.origin.y);
            var norm = v.norm();

            var theta = (v.x * this.edge.direction.x + v.y * this.edge.direction.y) / norm;
            var prog = norm * Math.abs(theta);
            // returns length of projection on edge
            return prog;
        }

    };
    return Ant;
}


},{"./initializePoints.js":"/Users/vallette/ANTS/AntColony/src/initializePoints.js","./mouse.js":"/Users/vallette/ANTS/AntColony/src/mouse.js","./utilities.js":"/Users/vallette/ANTS/AntColony/src/utilities.js","./vector.js":"/Users/vallette/ANTS/AntColony/src/vector.js"}],"/Users/vallette/ANTS/AntColony/src/createEdges.js":[function(require,module,exports){
'use strict'

var dt = require("delaunay-triangulate");

var range = require('./utilities.js').range;

var points = require('./initializePoints.js').points;
var textMesh = require('./initializePoints.js').textMesh;
//var textPointsId = require('./initializePoints.js').textPointsId;
var citySet = require('./initializePoints.js').citySet;
var nbRandomPoints = require('./initializePoints.js').nbRandomPoints;
var forcedEdges = require('./initializePoints.js').forcedEdges;

var Edge = require('./edge.js');

// triangulate
var cells = dt(points.map(function(p){
    return [p.x, p.y]
}))

var edges = [];
var permutations = [[0,1], [0,2], [1,2]];

// force the edges of the text to be edges of the graph
if (textMesh) {
    range(0, points.length - nbRandomPoints).forEach(function(id){
        var directLink = forcedEdges[id];
        var textEdge = Edge.create(points[id], points[directLink]);
        edges.push(textEdge);
        points[id].nexts.push(textEdge);
        //addToNextEdges(points[pt], undefined, [textEdge], nextEdges);
    })
}


cells.forEach(function(cell){
   
    for (var i = 0; i < 3; ++i){  // for each point.id listed in current cell
        var pt = points[cell[i]];
        //console.log("-------- id:" + pt);

        for (var j = 1; j < 3; ++j){ 

            var ptj = points[cell[( i + j ) % 3]]; // pick one of the other 2 points of the cell
            //console.log("--- other: " + ptj);
            var newEdge = undefined;

            // if pt already has nextEdges ...
            if (pt.nexts.length != 0) {
                //console.log("Trouvé 1");
                
                // ... get the points corresponding ...
                var tempPoints = pt.nexts.map(function(e){
                    return [e.pt1, e.pt2];
                }).reduce(function(a, b){
                     return a.concat(b);
                });

                //console.log(tempPoints);

                // ... and check if ptj already is part of the existing nextEdges. If not, add the edge.
                if (tempPoints.indexOf(ptj) == -1){
                    //console.log("Trouvé 2");
                    //newEdge = createEdge(points[pt], points[ptj]);
                    newEdge = Edge.create(pt, ptj);
                    edges.push(newEdge);
                    pt.nexts.push(newEdge);
                    //console.log(edges.length);
                }
            }
            else {
                //console.log("Trouvé 3");
                newEdge = Edge.create(pt, ptj);
                edges.push(newEdge);
                pt.nexts.push(newEdge);
                //console.log(edges.length);
            }

            // add also the edge to the edge's other point's nextEdges
            if (newEdge != undefined){
                ptj.nexts.push(newEdge);
                //addToNextEdges(points[ptj], cell, temp, nextEdges); 
            }         
        }

        // add the textEdges to nextEdges map
        //if (textMesh && textPointsId.indexOf(pt.id) != -1 && pt.id < (textPointsId.length - 1)) 
        // if (textMesh && citySet.has(pt))
        if (textMesh && (citySet.indexOf(pt.id) != -1))
        {
            // console.log('verif');
            var textEdge = Edge.create(pt, points[pt.id + 1]);
            //var temp2 = [];
            edges.push(textEdge);
            pt.nexts.push(textEdge);
            //temp2.push(textEdge);
            //addToNextEdges(points[pt+1], undefined, temp2, nextEdges);
        }

        //console.log("longueur: " + points[pt].nexts.length);
        //addToNextEdges(points[pt], cell, nexts, nextEdges);
    }
    //console.log(edges.length);
})

module.exports = edges;
},{"./edge.js":"/Users/vallette/ANTS/AntColony/src/edge.js","./initializePoints.js":"/Users/vallette/ANTS/AntColony/src/initializePoints.js","./utilities.js":"/Users/vallette/ANTS/AntColony/src/utilities.js","delaunay-triangulate":"/Users/vallette/ANTS/AntColony/node_modules/delaunay-triangulate/triangulate.js"}],"/Users/vallette/ANTS/AntColony/src/edge.js":[function(require,module,exports){
'use strict';

var sqrt = Math.sqrt;
var pow = Math.pow;
var abs = Math.abs;
var atan = Math.atan;

var Vector = require('./vector.js');


function Edge(ptA, ptB) {
    var distance = sqrt( pow(ptA.x - ptB.x, 2) + pow(ptA.y - ptB.y, 2) );

    // find line equation ax + by + c = 0
    var a = 1;
    var b = - (ptB.x - ptA.x) / (ptB.y - ptA.y);

    // orientate vector (a,b)
    if (b < 0){
        b = -b;
        a = -a;
    }

    // normalize vector (a,b)
    var n = new Vector(a, b);
    n.normalize();

    var c = - (a * ptA.x + b * ptA.y);

    // // calculate vector director
    var v = new Vector(ptB.x - ptA.x, ptB.y - ptA.y);
    
    v.normalize();


    this.id = undefined;
    this.pt1 = ptA;
    this.pt2 = ptB;
    this.direction = v;
    this.orthDirection = n; 
    this.distance = distance;
    this.pheromon = 1/distance;
    this.line = {
        a: a,
        b: b,
        c: c,
    }
}


// static methods
Edge.create = function(ptA, ptB) {
    var edge = new Edge(ptA, ptB);
    return edge;
}


// methods
Edge.prototype = {

    getOtherPoint: function(point) {
        if (point == this.pt1)
            return this.pt2;
        else if (point == this.pt2)
            return this.pt1;
        else
            console.log("Error");
    },

    calculateDistance: function(x, y) {
        var a = this.line.a,
            b = this.line.b,
            c = this.line.c;
        return abs(a * x + b * y + c) / Math.sqrt(Math.pow(a,2) + Math.pow(b,2));
    },

    /*calculateDistance: function(point) {
        var a = this.line.a,
            b = this.line.b,
            c = this.line.c;
        return abs(a * point.x + b * point.y + c) / Math.sqrt(Math.pow(a,2) + Math.pow(b,2));
    }*/
}
module.exports = Edge;
},{"./vector.js":"/Users/vallette/ANTS/AntColony/src/vector.js"}],"/Users/vallette/ANTS/AntColony/src/initializeAnts.js":[function(require,module,exports){
'use strict'

var antFunction = require('./ant.js');

var nbAnts = 3000;

module.exports = function (container) {

	var Ant = antFunction(container);

	var population = new Array(nbAnts);
	var possibleStartPointsId = require('./initializePoints.js').possibleStartPointsId;

	for (var i = 0; i < nbAnts; i++) {
	    var newAnt = new Ant(Ant.generateRandStartPoint());
	    newAnt.setDirection();
	    population[i] = newAnt;
	}

	return population;

}
},{"./ant.js":"/Users/vallette/ANTS/AntColony/src/ant.js","./initializePoints.js":"/Users/vallette/ANTS/AntColony/src/initializePoints.js"}],"/Users/vallette/ANTS/AntColony/src/initializePoints.js":[function(require,module,exports){
'use strict'

var parse = require('parse-svg-path');

var range = require('./utilities.js').range;

var Point = require('./point.js');

var random = Math.random;

var nbRandomPoints = 300;
var nbStartPoints = 20;

var nbCity = 2;

var textMesh = true;

// Frame definition
var xInit = 0, yInit = 0;
var w = 1,
    h = 1;

var Achar = "c 4.6011,-11.71047 9.20835,-23.42006 13.8199,-35.12898 4.61156,-11.70892 9.22741,-23.41718 13.84573,-35.125 4.61831,-11.70782 9.23908,-23.41519 13.86046,-35.12233 4.62138,-11.70714 9.24336,-23.41406 13.86411,-35.12097 4.62074,-11.70691 9.24025,-23.41381 13.85667,-35.12092 4.61641,-11.70712 9.22974,-23.41444 13.83814,-35.12218 4.60839,-11.70775 9.21186,-23.41592 13.80854,-35.12474 4.59668,-11.70881 9.18658,-23.418277 13.76785,-35.128603 12.99923,-3.357855 24.30069,-6.821144 33.86893,-4.46411 9.56825,2.357035 17.40328,10.534393 23.46967,30.457833 4.17598,10.62114 8.36031,21.23882 12.54942,31.85453 4.1891,10.61571 8.38298,21.22943 12.57807,31.84266 4.19509,10.61323 8.3914,21.22595 12.58535,31.83965 4.19395,10.61369 8.38555,21.22836 12.57123,31.84549 4.18568,10.61712 8.36544,21.2367 12.53572,31.86021 4.17028,10.6235 8.33108,21.25093 12.47883,31.88377 4.14775,10.63283 8.28245,21.27107 12.40053,31.9162 4.11808,10.64512 8.21955,21.29712 12.30085,31.95749 -5.90896,6.95561 -24.61617,1.11298 -35.59372,3 -16.75248,2.84859 -26.96421,-0.41416 -28.40628,-19 -3.17726,-8.42157 -6.35606,-16.87924 -9.47997,-25.34854 -3.12391,-8.46929 -6.19293,-16.9502 -9.15062,-25.41826 -16.46723,-0.80053 -35.17768,-2.74281 -53.33727,-3.11439 -18.1596,-0.37158 -35.76834,0.82753 -50.03214,6.30976 -4.15679,10.93124 -8.13806,21.9269 -12.15785,32.90546 -4.0198,10.97855 -8.07813,21.94001 -12.38903,32.80282 -9.52758,0.82747 -19.10457,0.96423 -28.69282,0.93363 -9.58824,-0.0306 -19.18773,-0.22855 -28.7603,-0.0705 l 0,-2.19791 z m 174,-109 c -3.17462,-9.51136 -6.44835,-18.99174 -9.7543,-28.46224 -6.81559,-19.51412 -6.24202,-25.37946 -12.3174,-43.60788 -3.16722,-9.51543 -13.60257,-32.32866 -16.49973,-41.92988 -6.08989,12.88576 -11.132,26.59272 -15.93415,40.47476 -4.80215,13.88203 -9.36435,27.93915 -14.49442,41.52524 -2.39425,12.05801 -22.8815,40.11116 2.18745,34 11.05256,-0.44486 22.52939,0.11061 33.85624,0.24955 11.32684,0.13895 22.50369,-0.13863 32.95631,-2.24955 z ";
var Nchar = "c 0,0 0,-23.83334 0,-35.75 0,-11.91667 0,-23.83334 0,-35.75 0,-11.91666 0,-23.83333 0,-35.75 0,-11.916667 0,-23.833335 0,-35.750003 14.64729,1.05313 31.03193,-2.08808 44.60679,1.55415 6.18076,7.610996 12.35438,15.231342 18.51973,22.861013 6.16535,7.62967 12.32244,15.26866 18.47014,22.91696 6.1477,7.64829 12.28602,15.30588 18.41383,22.97274 6.12781,7.66686 12.24512,15.343 18.35081,23.02838 6.10568,7.68538 12.19975,15.38001 18.28107,23.08386 6.08132,7.70384 12.1499,15.41691 18.20462,23.13918 6.05472,7.72226 12.09557,15.45372 18.12145,23.19435 6.02587,7.74063 12.03676,15.49043 18.03156,23.24937 0.4558,-15.4914 0.72655,-30.98725 0.88119,-46.48582 0.15464,-15.49858 0.19319,-30.99988 0.18458,-46.50218 -0.009,-15.5023 -0.0643,-31.0056 -0.0983,-46.50817 -0.034,-15.50258 -0.0461,-31.004432 0.0325,-46.503833 9.33333,0 18.66667,0 28,0 9.33333,0 18.66666,0 28,0 0,11.916667 0,23.833332 0,35.750003 0,11.91666 0,23.83333 0,35.75 0,11.91666 0,23.83333 0,35.75 0,11.91667 0,23.83333 0,35.75 0,11.91666 0,23.83333 0,35.75 0,11.91667 0,23.83334 0,35.75 0,11.91666 0,23.83333 0,35.75 0,11.91667 0,23.83334 0,35.75 -14.6441,-1.05348 -31.026,2.08847 -44.59764,-1.55416 -6.27247,-7.38037 -12.48436,-14.81692 -18.65123,-22.29507 -6.16688,-7.47815 -12.28873,-14.9979 -18.38111,-22.54465 -6.09238,-7.54675 -12.15528,-15.12051 -18.20425,-22.70667 -6.04896,-7.58617 -12.08399,-15.18476 -18.12063,-22.78116 -6.03664,-7.5964 -12.07489,-15.19062 -18.1303,-22.76807 -6.05541,-7.57745 -12.12797,-15.13813 -18.23323,-22.66744 -6.10526,-7.52932 -12.24322,-15.02727 -18.42942,-22.47926 -6.18619,-7.45199 -12.42063,-14.85803 -18.71886,-22.20352 -0.22791,15.16496 -0.36674,30.3308 -0.4489,45.49718 -0.0822,15.16637 -0.10766,30.33329 -0.10889,45.5004 -0.001,15.16712 0.0218,30.33443 0.0367,45.50162 0.0149,15.16718 0.0216,30.33422 -0.0122,45.5008 -9.33333,0 -18.66667,0 -28,0 -9.33333,0 -18.66666,0 -28,0 0,-11.91667 0,-23.83334 0,-35.75 0,-11.91667 0,-23.83334 0,-35.75 0,-11.91666 0,-23.83333 0,-35.75 0,-11.91667 0,-35.75 -10e-6,-35.75 z ";
var Tchar = "c 0,-9.91667 0,-19.83334 0,-29.75 0,-9.91667 0,-19.83334 0,-29.75 0,-9.91666 0,-19.83333 0,-29.75 0,-9.91666 -1.47394,-19.83333 -1.47394,-29.75 0,-15.33334 -29.19273,0 -44.52606,0 -15.33333,0 -30.66667,0 -46,0 0,-8 0,-16 0,-24.000002 0,-8 0,-16.000001 0,-24.000001 9.91666,0 19.83333,0 29.75,0 9.91667,0 19.83333,0 29.75,0 9.91666,0 19.83333,0 29.75,0 9.91667,0 19.83333,0 29.75,0 9.91666,0 19.83333,0 29.75,0 9.91667,0 19.83333,0 29.75,0 9.91666,0 19.83333,0 29.75,0 9.91667,0 19.83333,0 29.75,0 0,8 0,16.000001 0,24.000001 0,8.000002 0,16.000002 0,24.000002 -15,0 -30,0 -45,0 -15,0 -44.88809,-13.52565 -45,1.47394 -0.0735,9.85588 -0.10613,18.2399 -0.11249,28.09923 -0.006,9.85932 0.0135,19.72001 0.045,29.58133 0.0315,9.86132 0.0745,19.72329 0.1143,29.58517 0.0399,9.86188 0.0765,19.72367 0.0954,29.58467 0.0189,9.86099 0.0199,19.72118 -0.0117,29.57984 -0.0315,9.85866 -0.0956,19.7158 -0.20688,29.57069 -0.11131,9.85488 -0.26985,19.70752 -0.49033,29.55719 -0.22047,9.84967 -0.50289,19.69636 -0.86193,29.53937 -7.39624,-0.99964 -18.95658,0.96726 -29.70912,1.82971 -10.75253,0.86245 -24.48556,2.02609 -24.86231,-4.79696 -0.66223,-11.99314 -0.66223,-21.54349 -0.49667,-30.48314 0.16556,-8.93965 0.49667,-17.2686 0.49667,-26.81895 0,-9.55035 0,-19.1007 0,-28.65105 0,-9.55035 0,-19.10069 0,-28.65104 z";
var Schar = "c -16.6587,-2.34387 -33.35995,-6.23515 -49.25137,-12.08831 -15.89142,-5.85316 -30.97305,-13.6682 -44.3926,-23.85957 3.73715,-8.02639 7.96065,-15.79371 12.19101,-23.55228 4.23037,-7.75858 8.46759,-15.50842 12.23219,-23.49984 8.82902,6.69162 18.48235,12.75663 28.67412,17.86344 10.19178,5.10681 20.92205,9.25542 31.90485,12.11421 10.9829,2.8588 22.2183,4.4278 33.4206,4.37539 11.2022,-0.0524 22.3713,-1.72622 33.2212,-5.35304 14.2804,-6.92169 17.3033,-19.64436 13.4946,-31.15028 -3.8087,-11.50591 -14.4489,-21.79507 -27.4946,-23.84972 -10.162,-4.07686 -21.1052,-7.14863 -32.1975,-10.11097 -11.0924,-2.96235 -22.334,-5.81526 -33.0931,-9.45439 -10.7591,-3.63914 -21.03562,-8.06449 -30.19779,-14.17171 -9.16216,-6.10723 -17.20996,-13.89632 -23.51158,-24.26293 -5.1794,-10.97058 -7.48544,-22.77583 -7.30512,-34.52966 0.18031,-11.75382 2.84698,-23.4562 7.61299,-34.22103 4.76601,-10.76483 11.63137,-20.5921 20.20905,-28.595692 8.57769,-8.003592 18.86771,-14.183506 30.48305,-17.653621 12.1821,-4.024145 24.8777,-6.357009 37.6988,-7.121901 12.8211,-0.764892 25.7677,0.03819 38.4517,2.285933 12.6841,2.247744 25.1055,5.940153 36.8765,10.953918 11.7709,5.013764 22.8912,11.348885 32.973,18.882053 -4.4745,7.00368 -8.2581,14.70708 -12.1886,22.23499 -3.9306,7.5279 -8.008,14.88031 -13.0701,21.18198 -15.0728,-11.15158 -33.5009,-20.54491 -52.6114,-24.91726 -19.1105,-4.37235 -38.9034,-3.72372 -56.7058,5.20861 -10.0321,7.02789 -13.3348,18.84695 -11.1561,29.53596 2.1787,10.68902 9.8387,20.24799 21.732,22.75572 10.1335,4.43307 21.0241,7.69407 32.1101,10.69023 11.0861,2.99616 22.3676,5.72747 33.283,9.10117 10.9155,3.3737 21.4648,7.38978 31.0865,12.95547 9.6217,5.56569 18.3157,12.68099 25.5204,22.25313 6.6301,9.71635 10.5587,20.88023 12.0064,32.41231 1.4476,11.53207 0.4142,23.43234 -2.8797,34.62146 -3.2939,11.18912 -8.8484,21.66709 -16.443,30.35457 -7.5946,8.68749 -17.2293,15.58449 -28.6837,19.61166 -13.0091,5.92928 -27.0197,8.63869 -41.2728,9.63608 -14.253,0.99738 -28.7484,0.28274 -42.7272,-0.63608 z ";

var svgString = "m 0,0 " + Achar +"m 126,-31 " + Nchar + "m 376,24 "+ Tchar +"m 250,120 "+ Schar;

 
function svgToPoints(svgString) {
    var points = [];
    //var edges = [];

    // var edges = new Map(); // pointId -> pointId on border
    var edges = Object.create(null);

    var beginingPath;

    var X = 0;
    var Y = 0;
    var nbPoints = 0;
    var prevPoint;

    var commands = parse(svgString)
    for (var i=0; i<commands.length; i++){
        var command = commands[i];
        switch (command[0]) {
            case "m":
                X += command[1];
                Y += command[2];
                prevPoint = undefined;
                beginingPath = nbPoints;
                break;
            case "M":
                X = command[1];
                Y = command[2];
                prevPoint = undefined;
                beginingPath = nbPoints;
                break;  
            case "c":
                X += command[5];
                Y += command[6];
                points.push({id:nbPoints, x:X, y:Y});
                // nbPoints++;
                // if (prevPoint) {
                //     edges.push([prevPoint, nbPoints]);
                //     prevPoint = nbPoints;
                // }
                if (prevPoint != undefined) {
                    // edges.set(prevPoint, nbPoints);
                    edges[prevPoint] = nbPoints;
                }
                prevPoint = nbPoints;
                nbPoints++;
                break;
            case "z":
                // edges.set(prevPoint, beginingPath);
                edges[prevPoint] = nbPoints;
                beginingPath = undefined;
                prevPoint = undefined;
                break;    
        }
    }
    return {points : points, edges : edges};
}

// initialize points
var points = [];
var forcedEdges;
var citySet;

if (textMesh){

    var myText = svgToPoints(svgString);
    points = myText.points;
    forcedEdges = myText.edges;
    // citySet = new Set(range(0, points.length));
    citySet = range(0, points.length);

    var scaleX = 0.5;
    var scaleY = 0.2;
    var deltaX = 0.25;
    var deltaY = 0.3;

    // scale points to [0,1] + delta
    var maxX = Math.max.apply(Math, points.map(function(p){return p.x}));
    var minX = Math.min.apply(Math, points.map(function(p){return p.x}));
    var maxY = Math.max.apply(Math, points.map(function(p){return p.y}));
    var minY = Math.min.apply(Math, points.map(function(p){return p.y}));
    points = points.map(function(p){
        var x = scaleX * (p.x-minX)/(maxX-minX) + deltaX;
        var y = scaleY * (p.y-minY)/(maxY-minY) + deltaY;
        var newPoint = new Point(x, y);
        newPoint.id = p.id;

        return newPoint;
    });

    // only add random points
    var nbPoints = points.length;
    for(var i=0; i<nbRandomPoints; ++i) {

        //var x = random() * w + xInit;
        //var y = random() * h + yInit;

        var x = random();
        var y = random();

        var newPoint = new Point(x, y);
        newPoint.id = nbPoints;

        points.push(newPoint);

        nbPoints++;
    }

} else {
    //add random points

    var nbPoints = 0;
    for(var i=0; i<nbRandomPoints; ++i) {

        var x = random();
        var y = random();

        var newPoint = new Point(x, y);
        newPoint.id = nbPoints;

        points.push(newPoint);
        
        nbPoints++;
    }

    // citySet = new Set(range(0, nbCity));
    citySet = range(0, nbCity);
    console.log(citySet);
}


// initialize start points
var possibleStartPointsId = [];

for (var i = 0; i < nbStartPoints; i++){
    possibleStartPointsId.push(Math.floor(nbRandomPoints * random()));
}



module.exports = {
    textMesh: textMesh,
    points: points,
    citySet: citySet,
    possibleStartPointsId: possibleStartPointsId,
    nbRandomPoints: nbRandomPoints,
    forcedEdges: forcedEdges
}
},{"./point.js":"/Users/vallette/ANTS/AntColony/src/point.js","./utilities.js":"/Users/vallette/ANTS/AntColony/src/utilities.js","parse-svg-path":"/Users/vallette/ANTS/AntColony/node_modules/parse-svg-path/index.js"}],"/Users/vallette/ANTS/AntColony/src/mouse.js":[function(require,module,exports){
'use strict'

module.exports = function (container){

	var mouse = {
	    x: 0,
	    y: 0
	};

	container.addEventListener( 'mousemove', function(e){
	    var rect = container.getBoundingClientRect();
	    mouse.x = e.clientX / rect.width;
	    mouse.y = e.clientY / rect.height;
	});

	return mouse;

};

},{}],"/Users/vallette/ANTS/AntColony/src/point.js":[function(require,module,exports){
'use strict'

function Point(x, y) {
    this.id = undefined;                
    this.x = x;
    this.y = y;
    this.nexts = [];
}

module.exports = Point;
},{}],"/Users/vallette/ANTS/AntColony/src/rendering.js":[function(require,module,exports){
'use strict'

var random = Math.random;

module.exports = function(container){
    
    if(!container)
        throw new TypeError('Missing container');

    var points = require('./initializePoints.js').points;
    var citySet = require('./initializePoints.js').citySet;

    var edges = require('./createEdges.js');

    var population = require('./initializeAnts')(container);
    var nbAnts = population.length;
        
    var canvas = document.createElement("canvas");
    var rect = container.getBoundingClientRect();
    canvas.width = rect.width;
    canvas.height = rect.height;
    container.appendChild(canvas);
    
    var context = canvas.getContext("2d");
    

    function tick() {
        var w = canvas.width;
        var h = canvas.height;
        var mouse = [lastMouseMoveEvent.clientX/w, lastMouseMoveEvent.clientY/h];
        context.setTransform(w, 0, 0, h, 0, 0);
        context.fillStyle = "#fff";
        context.fillRect(0,0,w,h);

        // edges
        // context.strokeStyle = "#000";
        // for(var i=0; i < edges.length; ++i) {
        //     context.lineWidth = 0.0001;
        //     var edge = edges[i];
        //     // if (edge.pheromon != 0){
        //     //     context.lineWidth = Math.min(0.00001 * edge.pheromon, 0.01);
        //     // } else {
        //     //     context.lineWidth = 0.00001;
        //     // }
        //     context.beginPath();
        //     context.moveTo(points[edge.pt1.id].x, points[edge.pt1.id].y);
        //     context.lineTo(points[edge.pt2.id].x, points[edge.pt2.id].y);
        //     context.stroke();
        // }

        // // vertices
        // for(var i=0; i<points.length; ++i) {
        //     context.beginPath()
        //     var point = points[i];
        //     if (citySet.has(point.id)) {
        //         context.fillStyle = "#0101DF";
        //         context.arc(point.x, point.y, 0.006, 0, 2*Math.PI);
        //     }
        //     else {
        //         context.fillStyle = "#000";
        //         context.arc(points[i].x, points[i].y, 0.003, 0, 2*Math.PI);
        //     }
        //     context.closePath();
        //     context.fill();
        // }

        // move ants
        for (i = 0; i < nbAnts; i++) {
            population[i].transit();
        }

        // pheromon evaporation
        for (i = 0; i < edges.length; i++) {
            if(edges[i].pheromon > 0){
                edges[i].pheromon -= 0.0001;
            }
        }

        // ants
        for(var i=0; i<population.length; ++i) {
            context.beginPath()
            var x = population[i].x + 0.005*random();
            var y = population[i].y + 0.005*random();
            // var x = population[i].x;
            // var y = population[i].y;
            context.fillStyle = "black"
            context.fillRect(x, y, 0.0012, 0.0012);
            context.closePath();
            context.fill();
        }

    };
    
    var lastMouseMoveEvent = {
        clientX: 0,
        clientY: 0
    };
    
    container.addEventListener('mousemove', function(e){
        lastMouseMoveEvent = e;
    });
    
    var paused = false;
    
    function togglePlayPause(){
        paused = !paused;
        if(!paused)
            animate();
    }
    
    container.addEventListener('click', togglePlayPause);
    
    function animate(){
        tick();
        
        if(!paused)
            requestAnimationFrame(animate);
    };
    animate();
    
    return {
        togglePlayPause: togglePlayPause,
        // should be a getter/setter, but IE8
        getAntCount: function(){
            throw 'TODO';
        },
        setAntCount: function(){
            throw 'TODO';
        }
    }
}

},{"./createEdges.js":"/Users/vallette/ANTS/AntColony/src/createEdges.js","./initializeAnts":"/Users/vallette/ANTS/AntColony/src/initializeAnts.js","./initializePoints.js":"/Users/vallette/ANTS/AntColony/src/initializePoints.js"}],"/Users/vallette/ANTS/AntColony/src/runAnimation.js":[function(require,module,exports){
'use strict';

require('../index.js')(document.body);
},{"../index.js":"/Users/vallette/ANTS/AntColony/index.js"}],"/Users/vallette/ANTS/AntColony/src/utilities.js":[function(require,module,exports){
'use strict'

var sqrt = Math.sqrt;
var pow = Math.pow;

function sign(x) {
	return x ? x < 0 ? -1 : 1 : 0;
}

function range(start, count) {
    return Array.apply(0, Array(count)).map(function (element, index) {
    	return index + start
    });
}

function distance(a, b){
	return sqrt(pow(a.x - b.x, 2) + pow(a.y - b.y, 2));
}

function norm(v){
	return sqrt(pow(v.x, 2) + pow(v.y, 2));
}

module.exports = {
	sign: sign,
	range: range,
	distance: distance,
	norm: norm
}
},{}],"/Users/vallette/ANTS/AntColony/src/vector.js":[function(require,module,exports){
'use strict'

function Vector(x, y) {
    this.x = x;                
    this.y = y;
}

Vector.prototype.norm = function(){
	return Math.sqrt(this.x * this.x + this.y * this.y);
}

Vector.prototype.normalize = function(){
	var norm = this.norm();
	this.x = this.x / norm;
	this.y = this.y / norm;
}



module.exports = Vector;
},{}]},{},["/Users/vallette/ANTS/AntColony/src/runAnimation.js"])
//# sourceMappingURL=data:application/json;base64,eyJ2ZXJzaW9uIjozLCJzb3VyY2VzIjpbIi91c3IvbG9jYWwvbGliL25vZGVfbW9kdWxlcy93YXRjaGlmeS9ub2RlX21vZHVsZXMvYnJvd3NlcmlmeS9ub2RlX21vZHVsZXMvYnJvd3Nlci1wYWNrL19wcmVsdWRlLmpzIiwiL1VzZXJzL3ZhbGxldHRlL0FOVFMvQW50Q29sb255L2luZGV4LmpzIiwiL1VzZXJzL3ZhbGxldHRlL0FOVFMvQW50Q29sb255L25vZGVfbW9kdWxlcy9kZWxhdW5heS10cmlhbmd1bGF0ZS9ub2RlX21vZHVsZXMvaW5jcmVtZW50YWwtY29udmV4LWh1bGwvaWNoLmpzIiwiL1VzZXJzL3ZhbGxldHRlL0FOVFMvQW50Q29sb255L25vZGVfbW9kdWxlcy9kZWxhdW5heS10cmlhbmd1bGF0ZS9ub2RlX21vZHVsZXMvaW5jcmVtZW50YWwtY29udmV4LWh1bGwvbm9kZV9tb2R1bGVzL3JvYnVzdC1vcmllbnRhdGlvbi9ub2RlX21vZHVsZXMvcm9idXN0LXNjYWxlL25vZGVfbW9kdWxlcy90d28tc3VtL3R3by1zdW0uanMiLCIvVXNlcnMvdmFsbGV0dGUvQU5UUy9BbnRDb2xvbnkvbm9kZV9tb2R1bGVzL2RlbGF1bmF5LXRyaWFuZ3VsYXRlL25vZGVfbW9kdWxlcy9pbmNyZW1lbnRhbC1jb252ZXgtaHVsbC9ub2RlX21vZHVsZXMvcm9idXN0LW9yaWVudGF0aW9uL25vZGVfbW9kdWxlcy9yb2J1c3Qtc2NhbGUvcm9idXN0LXNjYWxlLmpzIiwiL1VzZXJzL3ZhbGxldHRlL0FOVFMvQW50Q29sb255L25vZGVfbW9kdWxlcy9kZWxhdW5heS10cmlhbmd1bGF0ZS9ub2RlX21vZHVsZXMvaW5jcmVtZW50YWwtY29udmV4LWh1bGwvbm9kZV9tb2R1bGVzL3JvYnVzdC1vcmllbnRhdGlvbi9ub2RlX21vZHVsZXMvcm9idXN0LXN1YnRyYWN0L3JvYnVzdC1kaWZmLmpzIiwiL1VzZXJzL3ZhbGxldHRlL0FOVFMvQW50Q29sb255L25vZGVfbW9kdWxlcy9kZWxhdW5heS10cmlhbmd1bGF0ZS9ub2RlX21vZHVsZXMvaW5jcmVtZW50YWwtY29udmV4LWh1bGwvbm9kZV9tb2R1bGVzL3JvYnVzdC1vcmllbnRhdGlvbi9ub2RlX21vZHVsZXMvcm9idXN0LXN1bS9yb2J1c3Qtc3VtLmpzIiwiL1VzZXJzL3ZhbGxldHRlL0FOVFMvQW50Q29sb255L25vZGVfbW9kdWxlcy9kZWxhdW5heS10cmlhbmd1bGF0ZS9ub2RlX21vZHVsZXMvaW5jcmVtZW50YWwtY29udmV4LWh1bGwvbm9kZV9tb2R1bGVzL3JvYnVzdC1vcmllbnRhdGlvbi9ub2RlX21vZHVsZXMvdHdvLXByb2R1Y3QvdHdvLXByb2R1Y3QuanMiLCIvVXNlcnMvdmFsbGV0dGUvQU5UUy9BbnRDb2xvbnkvbm9kZV9tb2R1bGVzL2RlbGF1bmF5LXRyaWFuZ3VsYXRlL25vZGVfbW9kdWxlcy9pbmNyZW1lbnRhbC1jb252ZXgtaHVsbC9ub2RlX21vZHVsZXMvcm9idXN0LW9yaWVudGF0aW9uL29yaWVudGF0aW9uLmpzIiwiL1VzZXJzL3ZhbGxldHRlL0FOVFMvQW50Q29sb255L25vZGVfbW9kdWxlcy9kZWxhdW5heS10cmlhbmd1bGF0ZS9ub2RlX21vZHVsZXMvaW5jcmVtZW50YWwtY29udmV4LWh1bGwvbm9kZV9tb2R1bGVzL3NpbXBsaWNpYWwtY29tcGxleC9ub2RlX21vZHVsZXMvYml0LXR3aWRkbGUvdHdpZGRsZS5qcyIsIi9Vc2Vycy92YWxsZXR0ZS9BTlRTL0FudENvbG9ueS9ub2RlX21vZHVsZXMvZGVsYXVuYXktdHJpYW5ndWxhdGUvbm9kZV9tb2R1bGVzL2luY3JlbWVudGFsLWNvbnZleC1odWxsL25vZGVfbW9kdWxlcy9zaW1wbGljaWFsLWNvbXBsZXgvbm9kZV9tb2R1bGVzL3VuaW9uLWZpbmQvaW5kZXguanMiLCIvVXNlcnMvdmFsbGV0dGUvQU5UUy9BbnRDb2xvbnkvbm9kZV9tb2R1bGVzL2RlbGF1bmF5LXRyaWFuZ3VsYXRlL25vZGVfbW9kdWxlcy9pbmNyZW1lbnRhbC1jb252ZXgtaHVsbC9ub2RlX21vZHVsZXMvc2ltcGxpY2lhbC1jb21wbGV4L3RvcG9sb2d5LmpzIiwiL1VzZXJzL3ZhbGxldHRlL0FOVFMvQW50Q29sb255L25vZGVfbW9kdWxlcy9kZWxhdW5heS10cmlhbmd1bGF0ZS9ub2RlX21vZHVsZXMvdW5pcS91bmlxLmpzIiwiL1VzZXJzL3ZhbGxldHRlL0FOVFMvQW50Q29sb255L25vZGVfbW9kdWxlcy9kZWxhdW5heS10cmlhbmd1bGF0ZS90cmlhbmd1bGF0ZS5qcyIsIi9Vc2Vycy92YWxsZXR0ZS9BTlRTL0FudENvbG9ueS9ub2RlX21vZHVsZXMvcGFyc2Utc3ZnLXBhdGgvaW5kZXguanMiLCIvVXNlcnMvdmFsbGV0dGUvQU5UUy9BbnRDb2xvbnkvc3JjL2FudC5qcyIsIi9Vc2Vycy92YWxsZXR0ZS9BTlRTL0FudENvbG9ueS9zcmMvY3JlYXRlRWRnZXMuanMiLCIvVXNlcnMvdmFsbGV0dGUvQU5UUy9BbnRDb2xvbnkvc3JjL2VkZ2UuanMiLCIvVXNlcnMvdmFsbGV0dGUvQU5UUy9BbnRDb2xvbnkvc3JjL2luaXRpYWxpemVBbnRzLmpzIiwiL1VzZXJzL3ZhbGxldHRlL0FOVFMvQW50Q29sb255L3NyYy9pbml0aWFsaXplUG9pbnRzLmpzIiwiL1VzZXJzL3ZhbGxldHRlL0FOVFMvQW50Q29sb255L3NyYy9tb3VzZS5qcyIsIi9Vc2Vycy92YWxsZXR0ZS9BTlRTL0FudENvbG9ueS9zcmMvcG9pbnQuanMiLCIvVXNlcnMvdmFsbGV0dGUvQU5UUy9BbnRDb2xvbnkvc3JjL3JlbmRlcmluZy5qcyIsIi9Vc2Vycy92YWxsZXR0ZS9BTlRTL0FudENvbG9ueS9zcmMvcnVuQW5pbWF0aW9uLmpzIiwiL1VzZXJzL3ZhbGxldHRlL0FOVFMvQW50Q29sb255L3NyYy91dGlsaXRpZXMuanMiLCIvVXNlcnMvdmFsbGV0dGUvQU5UUy9BbnRDb2xvbnkvc3JjL3ZlY3Rvci5qcyJdLCJuYW1lcyI6W10sIm1hcHBpbmdzIjoiQUFBQTtBQ0FBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTs7QUNSQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBOztBQzdiQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBOztBQ2hCQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBOztBQ2pEQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7O0FDM0pBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTs7QUMzSkE7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBOztBQ2hDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTs7QUM3TEE7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7O0FDNU1BO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBOztBQ3pEQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTs7QUN0VkE7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7O0FDekRBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTs7QUM5SkE7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTs7QUN2REE7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTs7QUN2TUE7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7O0FDekdBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTs7QUNuRkE7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7O0FDckJBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTs7QUNoTEE7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7O0FDbEJBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBOztBQ1RBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTs7QUNuSUE7QUFDQTtBQUNBOztBQ0ZBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7O0FDNUJBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0EiLCJmaWxlIjoiZ2VuZXJhdGVkLmpzIiwic291cmNlUm9vdCI6IiIsInNvdXJjZXNDb250ZW50IjpbIihmdW5jdGlvbiBlKHQsbixyKXtmdW5jdGlvbiBzKG8sdSl7aWYoIW5bb10pe2lmKCF0W29dKXt2YXIgYT10eXBlb2YgcmVxdWlyZT09XCJmdW5jdGlvblwiJiZyZXF1aXJlO2lmKCF1JiZhKXJldHVybiBhKG8sITApO2lmKGkpcmV0dXJuIGkobywhMCk7dmFyIGY9bmV3IEVycm9yKFwiQ2Fubm90IGZpbmQgbW9kdWxlICdcIitvK1wiJ1wiKTt0aHJvdyBmLmNvZGU9XCJNT0RVTEVfTk9UX0ZPVU5EXCIsZn12YXIgbD1uW29dPXtleHBvcnRzOnt9fTt0W29dWzBdLmNhbGwobC5leHBvcnRzLGZ1bmN0aW9uKGUpe3ZhciBuPXRbb11bMV1bZV07cmV0dXJuIHMobj9uOmUpfSxsLGwuZXhwb3J0cyxlLHQsbixyKX1yZXR1cm4gbltvXS5leHBvcnRzfXZhciBpPXR5cGVvZiByZXF1aXJlPT1cImZ1bmN0aW9uXCImJnJlcXVpcmU7Zm9yKHZhciBvPTA7bzxyLmxlbmd0aDtvKyspcyhyW29dKTtyZXR1cm4gc30pIiwiJ3VzZSBzdHJpY3QnO1xuXG52YXIgaW5pdFJlbmRlcmluZyA9IHJlcXVpcmUoJy4vc3JjL3JlbmRlcmluZy5qcycpO1xuXG5tb2R1bGUuZXhwb3J0cyA9IGZ1bmN0aW9uKGNvbnRhaW5lckVsZW1lbnQpe1xuICAgIGluaXRSZW5kZXJpbmcoY29udGFpbmVyRWxlbWVudCk7XG4gICAgdmFyIHBvaW50cyA9IHJlcXVpcmUoJy4vc3JjL2luaXRpYWxpemVQb2ludHMuanMnKTtcbiAgICB2YXIgZWRnZXMgPSByZXF1aXJlKCcuL3NyYy9jcmVhdGVFZGdlcy5qcycpO1xufTsiLCJcInVzZSBzdHJpY3RcIlxuXG4vL0hpZ2ggbGV2ZWwgaWRlYTpcbi8vIDEuIFVzZSBDbGFya3NvbidzIGluY3JlbWVudGFsIGNvbnN0cnVjdGlvbiB0byBmaW5kIGNvbnZleCBodWxsXG4vLyAyLiBQb2ludCBsb2NhdGlvbiBpbiB0cmlhbmd1bGF0aW9uIGJ5IGp1bXAgYW5kIHdhbGtcblxubW9kdWxlLmV4cG9ydHMgPSBpbmNyZW1lbnRhbENvbnZleEh1bGxcblxudmFyIG9yaWVudCA9IHJlcXVpcmUoXCJyb2J1c3Qtb3JpZW50YXRpb25cIilcbnZhciBjb21wYXJlQ2VsbCA9IHJlcXVpcmUoXCJzaW1wbGljaWFsLWNvbXBsZXhcIikuY29tcGFyZUNlbGxzXG5cbmZ1bmN0aW9uIGNvbXBhcmVJbnQoYSwgYikge1xuICByZXR1cm4gYSAtIGJcbn1cblxuZnVuY3Rpb24gU2ltcGxleCh2ZXJ0aWNlcywgYWRqYWNlbnQsIGJvdW5kYXJ5KSB7XG4gIHRoaXMudmVydGljZXMgPSB2ZXJ0aWNlc1xuICB0aGlzLmFkamFjZW50ID0gYWRqYWNlbnRcbiAgdGhpcy5ib3VuZGFyeSA9IGJvdW5kYXJ5XG4gIHRoaXMubGFzdFZpc2l0ZWQgPSAtMVxufVxuXG5TaW1wbGV4LnByb3RvdHlwZS5mbGlwID0gZnVuY3Rpb24oKSB7XG4gIHZhciB0ID0gdGhpcy52ZXJ0aWNlc1swXVxuICB0aGlzLnZlcnRpY2VzWzBdID0gdGhpcy52ZXJ0aWNlc1sxXVxuICB0aGlzLnZlcnRpY2VzWzFdID0gdFxuICB2YXIgdSA9IHRoaXMuYWRqYWNlbnRbMF1cbiAgdGhpcy5hZGphY2VudFswXSA9IHRoaXMuYWRqYWNlbnRbMV1cbiAgdGhpcy5hZGphY2VudFsxXSA9IHVcbn1cblxuZnVuY3Rpb24gR2x1ZUZhY2V0KHZlcnRpY2VzLCBjZWxsLCBpbmRleCkge1xuICB0aGlzLnZlcnRpY2VzID0gdmVydGljZXNcbiAgdGhpcy5jZWxsID0gY2VsbFxuICB0aGlzLmluZGV4ID0gaW5kZXhcbn1cblxuZnVuY3Rpb24gY29tcGFyZUdsdWUoYSwgYikge1xuICByZXR1cm4gY29tcGFyZUNlbGwoYS52ZXJ0aWNlcywgYi52ZXJ0aWNlcylcbn1cblxuZnVuY3Rpb24gYmFrZU9yaWVudChkKSB7XG4gIHZhciBjb2RlID0gW1wiZnVuY3Rpb24gb3JpZW50KCl7dmFyIHR1cGxlPXRoaXMudHVwbGU7cmV0dXJuIHRlc3QoXCJdXG4gIGZvcih2YXIgaT0wOyBpPD1kOyArK2kpIHtcbiAgICBpZihpID4gMCkge1xuICAgICAgY29kZS5wdXNoKFwiLFwiKVxuICAgIH1cbiAgICBjb2RlLnB1c2goXCJ0dXBsZVtcIiwgaSwgXCJdXCIpXG4gIH1cbiAgY29kZS5wdXNoKFwiKX1yZXR1cm4gb3JpZW50XCIpXG4gIHZhciBwcm9jID0gbmV3IEZ1bmN0aW9uKFwidGVzdFwiLCBjb2RlLmpvaW4oXCJcIikpXG4gIHZhciB0ZXN0ID0gb3JpZW50W2QrMV1cbiAgaWYoIXRlc3QpIHtcbiAgICB0ZXN0ID0gb3JpZW50XG4gIH1cbiAgcmV0dXJuIHByb2ModGVzdClcbn1cblxudmFyIEJBS0VEID0gW11cblxuZnVuY3Rpb24gVHJpYW5ndWxhdGlvbihkaW1lbnNpb24sIHZlcnRpY2VzLCBzaW1wbGljZXMpIHtcbiAgdGhpcy5kaW1lbnNpb24gPSBkaW1lbnNpb25cbiAgdGhpcy52ZXJ0aWNlcyA9IHZlcnRpY2VzXG4gIHRoaXMuc2ltcGxpY2VzID0gc2ltcGxpY2VzXG4gIHRoaXMuaW50ZXJpb3IgPSBzaW1wbGljZXMuZmlsdGVyKGZ1bmN0aW9uKGMpIHtcbiAgICByZXR1cm4gIWMuYm91bmRhcnlcbiAgfSlcblxuICB0aGlzLnR1cGxlID0gbmV3IEFycmF5KGRpbWVuc2lvbisxKVxuICBmb3IodmFyIGk9MDsgaTw9ZGltZW5zaW9uOyArK2kpIHtcbiAgICB0aGlzLnR1cGxlW2ldID0gdGhpcy52ZXJ0aWNlc1tpXVxuICB9XG5cbiAgdmFyIG8gPSBCQUtFRFtkaW1lbnNpb25dXG4gIGlmKCFvKSB7XG4gICAgbyA9IEJBS0VEW2RpbWVuc2lvbl0gPSBiYWtlT3JpZW50KGRpbWVuc2lvbilcbiAgfVxuICB0aGlzLm9yaWVudCA9IG9cbn1cblxudmFyIHByb3RvID0gVHJpYW5ndWxhdGlvbi5wcm90b3R5cGVcblxuLy9EZWdlbmVyYXRlIHNpdHVhdGlvbiB3aGVyZSB3ZSBhcmUgb24gYm91bmRhcnksIGJ1dCBjb3BsYW5hciB0byBmYWNlXG5wcm90by5oYW5kbGVCb3VuZGFyeURlZ2VuZXJhY3kgPSBmdW5jdGlvbihjZWxsLCBwb2ludCkge1xuICB2YXIgZCA9IHRoaXMuZGltZW5zaW9uXG4gIHZhciBuID0gdGhpcy52ZXJ0aWNlcy5sZW5ndGggLSAxXG4gIHZhciB0dXBsZSA9IHRoaXMudHVwbGVcbiAgdmFyIHZlcnRzID0gdGhpcy52ZXJ0aWNlc1xuXG4gIC8vRHVtYiBzb2x1dGlvbjogSnVzdCBkbyBkZnMgZnJvbSBib3VuZGFyeSBjZWxsIHVudGlsIHdlIGZpbmQgYW55IHBlYWssIG9yIHRlcm1pbmF0ZVxuICB2YXIgdG9WaXNpdCA9IFsgY2VsbCBdXG4gIGNlbGwubGFzdFZpc2l0ZWQgPSAtblxuICB3aGlsZSh0b1Zpc2l0Lmxlbmd0aCA+IDApIHtcbiAgICBjZWxsID0gdG9WaXNpdC5wb3AoKVxuICAgIHZhciBjZWxsVmVydHMgPSBjZWxsLnZlcnRpY2VzXG4gICAgdmFyIGNlbGxBZGogPSBjZWxsLmFkamFjZW50XG4gICAgZm9yKHZhciBpPTA7IGk8PWQ7ICsraSkge1xuICAgICAgdmFyIG5laWdoYm9yID0gY2VsbEFkaltpXVxuICAgICAgaWYoIW5laWdoYm9yLmJvdW5kYXJ5IHx8IG5laWdoYm9yLmxhc3RWaXNpdGVkIDw9IC1uKSB7XG4gICAgICAgIGNvbnRpbnVlXG4gICAgICB9XG4gICAgICB2YXIgbnYgPSBuZWlnaGJvci52ZXJ0aWNlc1xuICAgICAgZm9yKHZhciBqPTA7IGo8PWQ7ICsraikge1xuICAgICAgICB2YXIgdnYgPSBudltqXVxuICAgICAgICBpZih2diA8IDApIHtcbiAgICAgICAgICB0dXBsZVtqXSA9IHBvaW50XG4gICAgICAgIH0gZWxzZSB7XG4gICAgICAgICAgdHVwbGVbal0gPSB2ZXJ0c1t2dl1cbiAgICAgICAgfVxuICAgICAgfVxuICAgICAgdmFyIG8gPSB0aGlzLm9yaWVudCgpXG4gICAgICBpZihvID4gMCkge1xuICAgICAgICByZXR1cm4gbmVpZ2hib3JcbiAgICAgIH1cbiAgICAgIG5laWdoYm9yLmxhc3RWaXNpdGVkID0gLW5cbiAgICAgIGlmKG8gPT09IDApIHtcbiAgICAgICAgdG9WaXNpdC5wdXNoKG5laWdoYm9yKVxuICAgICAgfVxuICAgIH1cbiAgfVxuICByZXR1cm4gbnVsbFxufVxuXG5wcm90by53YWxrID0gZnVuY3Rpb24ocG9pbnQsIHJhbmRvbSkge1xuICAvL0FsaWFzIGxvY2FsIHByb3BlcnRpZXNcbiAgdmFyIG4gPSB0aGlzLnZlcnRpY2VzLmxlbmd0aCAtIDFcbiAgdmFyIGQgPSB0aGlzLmRpbWVuc2lvblxuICB2YXIgdmVydHMgPSB0aGlzLnZlcnRpY2VzXG4gIHZhciB0dXBsZSA9IHRoaXMudHVwbGVcblxuICAvL0NvbXB1dGUgaW5pdGlhbCBqdW1wIGNlbGxcbiAgdmFyIGluaXRJbmRleCA9IHJhbmRvbSA/ICh0aGlzLmludGVyaW9yLmxlbmd0aCAqIE1hdGgucmFuZG9tKCkpfDAgOiAodGhpcy5pbnRlcmlvci5sZW5ndGgtMSlcbiAgdmFyIGNlbGwgPSB0aGlzLmludGVyaW9yWyBpbml0SW5kZXggXVxuXG4gIC8vU3RhcnQgd2Fsa2luZ1xub3V0ZXJMb29wOlxuICB3aGlsZSghY2VsbC5ib3VuZGFyeSkge1xuICAgIHZhciBjZWxsVmVydHMgPSBjZWxsLnZlcnRpY2VzXG4gICAgdmFyIGNlbGxBZGogPSBjZWxsLmFkamFjZW50XG5cbiAgICBmb3IodmFyIGk9MDsgaTw9ZDsgKytpKSB7XG4gICAgICB0dXBsZVtpXSA9IHZlcnRzW2NlbGxWZXJ0c1tpXV1cbiAgICB9XG4gICAgY2VsbC5sYXN0VmlzaXRlZCA9IG5cblxuICAgIC8vRmluZCBmYXJ0aGVzdCBhZGphY2VudCBjZWxsXG4gICAgZm9yKHZhciBpPTA7IGk8PWQ7ICsraSkge1xuICAgICAgdmFyIG5laWdoYm9yID0gY2VsbEFkaltpXVxuICAgICAgaWYobmVpZ2hib3IubGFzdFZpc2l0ZWQgPj0gbikge1xuICAgICAgICBjb250aW51ZVxuICAgICAgfVxuICAgICAgdmFyIHByZXYgPSB0dXBsZVtpXVxuICAgICAgdHVwbGVbaV0gPSBwb2ludFxuICAgICAgdmFyIG8gPSB0aGlzLm9yaWVudCgpXG4gICAgICB0dXBsZVtpXSA9IHByZXZcbiAgICAgIGlmKG8gPCAwKSB7XG4gICAgICAgIGNlbGwgPSBuZWlnaGJvclxuICAgICAgICBjb250aW51ZSBvdXRlckxvb3BcbiAgICAgIH0gZWxzZSB7XG4gICAgICAgIGlmKCFuZWlnaGJvci5ib3VuZGFyeSkge1xuICAgICAgICAgIG5laWdoYm9yLmxhc3RWaXNpdGVkID0gblxuICAgICAgICB9IGVsc2Uge1xuICAgICAgICAgIG5laWdoYm9yLmxhc3RWaXNpdGVkID0gLW5cbiAgICAgICAgfVxuICAgICAgfVxuICAgIH1cbiAgICByZXR1cm5cbiAgfVxuXG4gIHJldHVybiBjZWxsXG59XG5cbnByb3RvLmFkZFBlYWtzID0gZnVuY3Rpb24ocG9pbnQsIGNlbGwpIHtcbiAgdmFyIG4gPSB0aGlzLnZlcnRpY2VzLmxlbmd0aCAtIDFcbiAgdmFyIGQgPSB0aGlzLmRpbWVuc2lvblxuICB2YXIgdmVydHMgPSB0aGlzLnZlcnRpY2VzXG4gIHZhciB0dXBsZSA9IHRoaXMudHVwbGVcbiAgdmFyIGludGVyaW9yID0gdGhpcy5pbnRlcmlvclxuICB2YXIgc2ltcGxpY2VzID0gdGhpcy5zaW1wbGljZXNcblxuICAvL1dhbGtpbmcgZmluaXNoZWQgYXQgYm91bmRhcnksIHRpbWUgdG8gYWRkIHBlYWtzXG4gIHZhciB0b3Zpc2l0ID0gWyBjZWxsIF1cblxuICAvL1N0cmV0Y2ggaW5pdGlhbCBib3VuZGFyeSBjZWxsIGludG8gYSBwZWFrXG4gIGNlbGwubGFzdFZpc2l0ZWQgPSBuXG4gIGNlbGwudmVydGljZXNbY2VsbC52ZXJ0aWNlcy5pbmRleE9mKC0xKV0gPSBuXG4gIGNlbGwuYm91bmRhcnkgPSBmYWxzZVxuICBpbnRlcmlvci5wdXNoKGNlbGwpXG5cbiAgLy9SZWNvcmQgYSBsaXN0IG9mIGFsbCBuZXcgYm91bmRhcmllcyBjcmVhdGVkIGJ5IGFkZGVkIHBlYWtzIHNvIHdlIGNhbiBnbHVlIHRoZW0gdG9nZXRoZXIgd2hlbiB3ZSBhcmUgYWxsIGRvbmVcbiAgdmFyIGdsdWVGYWNldHMgPSBbXVxuXG4gIC8vRG8gYSB0cmF2ZXJzYWwgb2YgdGhlIGJvdW5kYXJ5IHdhbGtpbmcgb3V0d2FyZCBmcm9tIHN0YXJ0aW5nIHBlYWtcbiAgd2hpbGUodG92aXNpdC5sZW5ndGggPiAwKSB7XG4gICAgLy9Qb3Agb2ZmIHBlYWsgYW5kIHdhbGsgb3ZlciBhZGphY2VudCBjZWxsc1xuICAgIHZhciBjZWxsID0gdG92aXNpdC5wb3AoKVxuICAgIHZhciBjZWxsVmVydHMgPSBjZWxsLnZlcnRpY2VzXG4gICAgdmFyIGNlbGxBZGogPSBjZWxsLmFkamFjZW50XG4gICAgdmFyIGluZGV4T2ZOID0gY2VsbFZlcnRzLmluZGV4T2YobilcbiAgICBpZihpbmRleE9mTiA8IDApIHtcbiAgICAgIGNvbnRpbnVlXG4gICAgfVxuXG4gICAgZm9yKHZhciBpPTA7IGk8PWQ7ICsraSkge1xuICAgICAgaWYoaSA9PT0gaW5kZXhPZk4pIHtcbiAgICAgICAgY29udGludWVcbiAgICAgIH1cblxuICAgICAgLy9Gb3IgZWFjaCBib3VuZGFyeSBuZWlnaGJvciBvZiB0aGUgY2VsbFxuICAgICAgdmFyIG5laWdoYm9yID0gY2VsbEFkaltpXVxuICAgICAgaWYoIW5laWdoYm9yLmJvdW5kYXJ5IHx8IG5laWdoYm9yLmxhc3RWaXNpdGVkID49IG4pIHtcbiAgICAgICAgY29udGludWVcbiAgICAgIH1cblxuICAgICAgdmFyIG52ID0gbmVpZ2hib3IudmVydGljZXNcblxuICAgICAgLy9UZXN0IGlmIG5laWdoYm9yIGlzIGEgcGVha1xuICAgICAgaWYobmVpZ2hib3IubGFzdFZpc2l0ZWQgIT09IC1uKSB7ICAgICAgXG4gICAgICAgIC8vQ29tcHV0ZSBvcmllbnRhdGlvbiBvZiBwIHJlbGF0aXZlIHRvIGVhY2ggYm91bmRhcnkgcGVha1xuICAgICAgICB2YXIgaW5kZXhPZk5lZzEgPSAwXG4gICAgICAgIGZvcih2YXIgaj0wOyBqPD1kOyArK2opIHtcbiAgICAgICAgICBpZihudltqXSA8IDApIHtcbiAgICAgICAgICAgIGluZGV4T2ZOZWcxID0galxuICAgICAgICAgICAgdHVwbGVbal0gPSBwb2ludFxuICAgICAgICAgIH0gZWxzZSB7XG4gICAgICAgICAgICB0dXBsZVtqXSA9IHZlcnRzW252W2pdXVxuICAgICAgICAgIH1cbiAgICAgICAgfVxuICAgICAgICB2YXIgbyA9IHRoaXMub3JpZW50KClcblxuICAgICAgICAvL1Rlc3QgaWYgbmVpZ2hib3IgY2VsbCBpcyBhbHNvIGEgcGVha1xuICAgICAgICBpZihvID4gMCkge1xuICAgICAgICAgIG52W2luZGV4T2ZOZWcxXSA9IG5cbiAgICAgICAgICBuZWlnaGJvci5ib3VuZGFyeSA9IGZhbHNlXG4gICAgICAgICAgaW50ZXJpb3IucHVzaChuZWlnaGJvcilcbiAgICAgICAgICB0b3Zpc2l0LnB1c2gobmVpZ2hib3IpXG4gICAgICAgICAgbmVpZ2hib3IubGFzdFZpc2l0ZWQgPSBuXG4gICAgICAgICAgY29udGludWVcbiAgICAgICAgfSBlbHNlIHtcbiAgICAgICAgICBuZWlnaGJvci5sYXN0VmlzaXRlZCA9IC1uXG4gICAgICAgIH1cbiAgICAgIH1cblxuICAgICAgdmFyIG5hID0gbmVpZ2hib3IuYWRqYWNlbnRcblxuICAgICAgLy9PdGhlcndpc2UsIHJlcGxhY2UgbmVpZ2hib3Igd2l0aCBuZXcgZmFjZVxuICAgICAgdmFyIHZ2ZXJ0cyA9IGNlbGxWZXJ0cy5zbGljZSgpXG4gICAgICB2YXIgdmFkaiA9IGNlbGxBZGouc2xpY2UoKVxuICAgICAgdmFyIG5jZWxsID0gbmV3IFNpbXBsZXgodnZlcnRzLCB2YWRqLCB0cnVlKVxuICAgICAgc2ltcGxpY2VzLnB1c2gobmNlbGwpXG5cbiAgICAgIC8vQ29ubmVjdCB0byBuZWlnaGJvclxuICAgICAgdmFyIG9wcG9zaXRlID0gbmEuaW5kZXhPZihjZWxsKVxuICAgICAgaWYob3Bwb3NpdGUgPCAwKSB7XG4gICAgICAgIGNvbnRpbnVlXG4gICAgICB9XG4gICAgICBuYVtvcHBvc2l0ZV0gPSBuY2VsbFxuICAgICAgdmFkaltpbmRleE9mTl0gPSBuZWlnaGJvclxuXG4gICAgICAvL0Nvbm5lY3QgdG8gY2VsbFxuICAgICAgdnZlcnRzW2ldID0gLTFcbiAgICAgIHZhZGpbaV0gPSBjZWxsXG4gICAgICBjZWxsQWRqW2ldID0gbmNlbGxcblxuICAgICAgLy9GbGlwIGZhY2V0XG4gICAgICBuY2VsbC5mbGlwKClcblxuICAgICAgLy9BZGQgdG8gZ2x1ZSBsaXN0XG4gICAgICBmb3IodmFyIGo9MDsgajw9ZDsgKytqKSB7XG4gICAgICAgIHZhciB1dSA9IHZ2ZXJ0c1tqXVxuICAgICAgICBpZih1dSA8IDAgfHwgdXUgPT09IG4pIHtcbiAgICAgICAgICBjb250aW51ZVxuICAgICAgICB9XG4gICAgICAgIHZhciBuZmFjZSA9IG5ldyBBcnJheShkLTEpXG4gICAgICAgIHZhciBucHRyID0gMFxuICAgICAgICBmb3IodmFyIGs9MDsgazw9ZDsgKytrKSB7XG4gICAgICAgICAgdmFyIHZ2ID0gdnZlcnRzW2tdXG4gICAgICAgICAgaWYodnYgPCAwIHx8IGsgPT09IGopIHtcbiAgICAgICAgICAgIGNvbnRpbnVlXG4gICAgICAgICAgfVxuICAgICAgICAgIG5mYWNlW25wdHIrK10gPSB2dlxuICAgICAgICB9XG4gICAgICAgIGdsdWVGYWNldHMucHVzaChuZXcgR2x1ZUZhY2V0KG5mYWNlLCBuY2VsbCwgaikpXG4gICAgICB9XG4gICAgfVxuICB9XG5cbiAgLy9HbHVlIGJvdW5kYXJ5IGZhY2V0cyB0b2dldGhlclxuICBnbHVlRmFjZXRzLnNvcnQoY29tcGFyZUdsdWUpXG5cbiAgZm9yKHZhciBpPTA7IGkrMTxnbHVlRmFjZXRzLmxlbmd0aDsgaSs9Mikge1xuICAgIHZhciBhID0gZ2x1ZUZhY2V0c1tpXVxuICAgIHZhciBiID0gZ2x1ZUZhY2V0c1tpKzFdXG4gICAgdmFyIGFpID0gYS5pbmRleFxuICAgIHZhciBiaSA9IGIuaW5kZXhcbiAgICBpZihhaSA8IDAgfHwgYmkgPCAwKSB7XG4gICAgICBjb250aW51ZVxuICAgIH1cbiAgICBhLmNlbGwuYWRqYWNlbnRbYS5pbmRleF0gPSBiLmNlbGxcbiAgICBiLmNlbGwuYWRqYWNlbnRbYi5pbmRleF0gPSBhLmNlbGxcbiAgfVxufVxuXG5wcm90by5pbnNlcnQgPSBmdW5jdGlvbihwb2ludCwgcmFuZG9tKSB7XG4gIC8vQWRkIHBvaW50XG4gIHZhciB2ZXJ0cyA9IHRoaXMudmVydGljZXNcbiAgdmVydHMucHVzaChwb2ludClcblxuICB2YXIgY2VsbCA9IHRoaXMud2Fsayhwb2ludCwgcmFuZG9tKVxuICBpZighY2VsbCkge1xuICAgIHJldHVyblxuICB9XG5cbiAgLy9BbGlhcyBsb2NhbCBwcm9wZXJ0aWVzXG4gIHZhciBkID0gdGhpcy5kaW1lbnNpb25cbiAgdmFyIHR1cGxlID0gdGhpcy50dXBsZVxuXG4gIC8vRGVnZW5lcmF0ZSBjYXNlOiBJZiBwb2ludCBpcyBjb3BsYW5hciB0byBjZWxsLCB0aGVuIHdhbGsgdW50aWwgd2UgZmluZCBhIG5vbi1kZWdlbmVyYXRlIGJvdW5kYXJ5XG4gIGZvcih2YXIgaT0wOyBpPD1kOyArK2kpIHtcbiAgICB2YXIgdnYgPSBjZWxsLnZlcnRpY2VzW2ldXG4gICAgaWYodnYgPCAwKSB7XG4gICAgICB0dXBsZVtpXSA9IHBvaW50XG4gICAgfSBlbHNlIHtcbiAgICAgIHR1cGxlW2ldID0gdmVydHNbdnZdXG4gICAgfVxuICB9XG4gIHZhciBvID0gdGhpcy5vcmllbnQodHVwbGUpXG4gIGlmKG8gPCAwKSB7XG4gICAgcmV0dXJuXG4gIH0gZWxzZSBpZihvID09PSAwKSB7XG4gICAgY2VsbCA9IHRoaXMuaGFuZGxlQm91bmRhcnlEZWdlbmVyYWN5KGNlbGwsIHBvaW50KVxuICAgIGlmKCFjZWxsKSB7XG4gICAgICByZXR1cm5cbiAgICB9XG4gIH1cblxuICAvL0FkZCBwZWFrc1xuICB0aGlzLmFkZFBlYWtzKHBvaW50LCBjZWxsKVxufVxuXG4vL0V4dHJhY3QgYWxsIGJvdW5kYXJ5IGNlbGxzXG5wcm90by5ib3VuZGFyeSA9IGZ1bmN0aW9uKCkge1xuICB2YXIgZCA9IHRoaXMuZGltZW5zaW9uXG4gIHZhciBib3VuZGFyeSA9IFtdXG4gIHZhciBjZWxscyA9IHRoaXMuc2ltcGxpY2VzXG4gIHZhciBuYyA9IGNlbGxzLmxlbmd0aFxuICBmb3IodmFyIGk9MDsgaTxuYzsgKytpKSB7XG4gICAgdmFyIGMgPSBjZWxsc1tpXVxuICAgIGlmKGMuYm91bmRhcnkpIHtcbiAgICAgIHZhciBiY2VsbCA9IG5ldyBBcnJheShkKVxuICAgICAgdmFyIGN2ID0gYy52ZXJ0aWNlc1xuICAgICAgdmFyIHB0ciA9IDBcbiAgICAgIHZhciBwYXJpdHkgPSAwXG4gICAgICBmb3IodmFyIGo9MDsgajw9ZDsgKytqKSB7XG4gICAgICAgIGlmKGN2W2pdID49IDApIHtcbiAgICAgICAgICBiY2VsbFtwdHIrK10gPSBjdltqXVxuICAgICAgICB9IGVsc2Uge1xuICAgICAgICAgIHBhcml0eSA9IGomMVxuICAgICAgICB9XG4gICAgICB9XG4gICAgICBpZihwYXJpdHkgPT09IChkJjEpKSB7XG4gICAgICAgIHZhciB0ID0gYmNlbGxbMF1cbiAgICAgICAgYmNlbGxbMF0gPSBiY2VsbFsxXVxuICAgICAgICBiY2VsbFsxXSA9IHRcbiAgICAgIH1cbiAgICAgIGJvdW5kYXJ5LnB1c2goYmNlbGwpXG4gICAgfVxuICB9XG4gIHJldHVybiBib3VuZGFyeVxufVxuXG5mdW5jdGlvbiBpbmNyZW1lbnRhbENvbnZleEh1bGwocG9pbnRzLCByYW5kb21TZWFyY2gpIHtcbiAgdmFyIG4gPSBwb2ludHMubGVuZ3RoXG4gIGlmKG4gPT09IDApIHtcbiAgICB0aHJvdyBuZXcgRXJyb3IoXCJNdXN0IGhhdmUgYXQgbGVhc3QgZCsxIHBvaW50c1wiKVxuICB9XG4gIHZhciBkID0gcG9pbnRzWzBdLmxlbmd0aFxuICBpZihuIDw9IGQpIHtcbiAgICB0aHJvdyBuZXcgRXJyb3IoXCJNdXN0IGlucHV0IGF0IGxlYXN0IGQrMSBwb2ludHNcIilcbiAgfVxuXG4gIC8vRklYTUU6IFRoaXMgY291bGQgYmUgZGVnZW5lcmF0ZSwgYnV0IG5lZWQgdG8gc2VsZWN0IGQrMSBub24tY29wbGFuYXIgcG9pbnRzIHRvIGJvb3RzdHJhcCBwcm9jZXNzXG4gIHZhciBpbml0aWFsU2ltcGxleCA9IHBvaW50cy5zbGljZSgwLCBkKzEpXG5cbiAgLy9NYWtlIHN1cmUgaW5pdGlhbCBzaW1wbGV4IGlzIHBvc2l0aXZlbHkgb3JpZW50ZWRcbiAgdmFyIG8gPSBvcmllbnQuYXBwbHkodm9pZCAwLCBpbml0aWFsU2ltcGxleClcbiAgaWYobyA9PT0gMCkge1xuICAgIHRocm93IG5ldyBFcnJvcihcIklucHV0IG5vdCBpbiBnZW5lcmFsIHBvc2l0aW9uXCIpXG4gIH1cbiAgdmFyIGluaXRpYWxDb29yZHMgPSBuZXcgQXJyYXkoZCsxKVxuICBmb3IodmFyIGk9MDsgaTw9ZDsgKytpKSB7XG4gICAgaW5pdGlhbENvb3Jkc1tpXSA9IGlcbiAgfVxuICBpZihvIDwgMCkge1xuICAgIGluaXRpYWxDb29yZHNbMF0gPSAxXG4gICAgaW5pdGlhbENvb3Jkc1sxXSA9IDBcbiAgfVxuXG4gIC8vQ3JlYXRlIGluaXRpYWwgdG9wb2xvZ2ljYWwgaW5kZXgsIGdsdWUgcG9pbnRlcnMgdG9nZXRoZXIgKGtpbmQgb2YgbWVzc3kpXG4gIHZhciBpbml0aWFsQ2VsbCA9IG5ldyBTaW1wbGV4KGluaXRpYWxDb29yZHMsIG5ldyBBcnJheShkKzEpLCBmYWxzZSlcbiAgdmFyIGJvdW5kYXJ5ID0gaW5pdGlhbENlbGwuYWRqYWNlbnRcbiAgdmFyIGxpc3QgPSBuZXcgQXJyYXkoZCsyKVxuICBmb3IodmFyIGk9MDsgaTw9ZDsgKytpKSB7XG4gICAgdmFyIHZlcnRzID0gaW5pdGlhbENvb3Jkcy5zbGljZSgpXG4gICAgZm9yKHZhciBqPTA7IGo8PWQ7ICsraikge1xuICAgICAgaWYoaiA9PT0gaSkge1xuICAgICAgICB2ZXJ0c1tqXSA9IC0xXG4gICAgICB9XG4gICAgfVxuICAgIHZhciB0ID0gdmVydHNbMF1cbiAgICB2ZXJ0c1swXSA9IHZlcnRzWzFdXG4gICAgdmVydHNbMV0gPSB0XG4gICAgdmFyIGNlbGwgPSBuZXcgU2ltcGxleCh2ZXJ0cywgbmV3IEFycmF5KGQrMSksIHRydWUpXG4gICAgYm91bmRhcnlbaV0gPSBjZWxsXG4gICAgbGlzdFtpXSA9IGNlbGxcbiAgfVxuICBsaXN0W2QrMV0gPSBpbml0aWFsQ2VsbFxuICBmb3IodmFyIGk9MDsgaTw9ZDsgKytpKSB7XG4gICAgdmFyIHZlcnRzID0gYm91bmRhcnlbaV0udmVydGljZXNcbiAgICB2YXIgYWRqID0gYm91bmRhcnlbaV0uYWRqYWNlbnRcbiAgICBmb3IodmFyIGo9MDsgajw9ZDsgKytqKSB7XG4gICAgICB2YXIgdiA9IHZlcnRzW2pdXG4gICAgICBpZih2IDwgMCkge1xuICAgICAgICBhZGpbal0gPSBpbml0aWFsQ2VsbFxuICAgICAgICBjb250aW51ZVxuICAgICAgfVxuICAgICAgZm9yKHZhciBrPTA7IGs8PWQ7ICsraykge1xuICAgICAgICBpZihib3VuZGFyeVtrXS52ZXJ0aWNlcy5pbmRleE9mKHYpIDwgMCkge1xuICAgICAgICAgIGFkaltqXSA9IGJvdW5kYXJ5W2tdXG4gICAgICAgIH1cbiAgICAgIH1cbiAgICB9XG4gIH1cblxuICAvL0luaXRpYWxpemUgdHJpYW5nbGVzXG4gIHZhciB0cmlhbmdsZXMgPSBuZXcgVHJpYW5ndWxhdGlvbihkLCBpbml0aWFsU2ltcGxleCwgbGlzdClcblxuICAvL0luc2VydCByZW1haW5pbmcgcG9pbnRzXG4gIHZhciB1c2VSYW5kb20gPSAhIXJhbmRvbVNlYXJjaFxuICBmb3IodmFyIGk9ZCsxOyBpPG47ICsraSkge1xuICAgIHRyaWFuZ2xlcy5pbnNlcnQocG9pbnRzW2ldLCB1c2VSYW5kb20pXG4gIH1cbiAgXG4gIC8vRXh0cmFjdCBib3VuZGFyeSBjZWxsc1xuICByZXR1cm4gdHJpYW5nbGVzLmJvdW5kYXJ5KClcbn0iLCJcInVzZSBzdHJpY3RcIlxuXG5tb2R1bGUuZXhwb3J0cyA9IGZhc3RUd29TdW1cblxuZnVuY3Rpb24gZmFzdFR3b1N1bShhLCBiLCByZXN1bHQpIHtcblx0dmFyIHggPSBhICsgYlxuXHR2YXIgYnYgPSB4IC0gYVxuXHR2YXIgYXYgPSB4IC0gYnZcblx0dmFyIGJyID0gYiAtIGJ2XG5cdHZhciBhciA9IGEgLSBhdlxuXHRpZihyZXN1bHQpIHtcblx0XHRyZXN1bHRbMF0gPSBhciArIGJyXG5cdFx0cmVzdWx0WzFdID0geFxuXHRcdHJldHVybiByZXN1bHRcblx0fVxuXHRyZXR1cm4gW2FyK2JyLCB4XVxufSIsIlwidXNlIHN0cmljdFwiXG5cbnZhciB0d29Qcm9kdWN0ID0gcmVxdWlyZShcInR3by1wcm9kdWN0XCIpXG52YXIgdHdvU3VtID0gcmVxdWlyZShcInR3by1zdW1cIilcblxubW9kdWxlLmV4cG9ydHMgPSBzY2FsZUxpbmVhckV4cGFuc2lvblxuXG5mdW5jdGlvbiBzY2FsZUxpbmVhckV4cGFuc2lvbihlLCBzY2FsZSkge1xuICB2YXIgbiA9IGUubGVuZ3RoXG4gIGlmKG4gPT09IDEpIHtcbiAgICB2YXIgdHMgPSB0d29Qcm9kdWN0KGVbMF0sIHNjYWxlKVxuICAgIGlmKHRzWzBdKSB7XG4gICAgICByZXR1cm4gdHNcbiAgICB9XG4gICAgcmV0dXJuIFsgdHNbMV0gXVxuICB9XG4gIHZhciBnID0gbmV3IEFycmF5KDIgKiBuKVxuICB2YXIgcSA9IFswLjEsIDAuMV1cbiAgdmFyIHQgPSBbMC4xLCAwLjFdXG4gIHZhciBjb3VudCA9IDBcbiAgdHdvUHJvZHVjdChlWzBdLCBzY2FsZSwgcSlcbiAgaWYocVswXSkge1xuICAgIGdbY291bnQrK10gPSBxWzBdXG4gIH1cbiAgZm9yKHZhciBpPTE7IGk8bjsgKytpKSB7XG4gICAgdHdvUHJvZHVjdChlW2ldLCBzY2FsZSwgdClcbiAgICB2YXIgcHEgPSBxWzFdXG4gICAgdHdvU3VtKHBxLCB0WzBdLCBxKVxuICAgIGlmKHFbMF0pIHtcbiAgICAgIGdbY291bnQrK10gPSBxWzBdXG4gICAgfVxuICAgIHZhciBhID0gdFsxXVxuICAgIHZhciBiID0gcVsxXVxuICAgIHZhciB4ID0gYSArIGJcbiAgICB2YXIgYnYgPSB4IC0gYVxuICAgIHZhciB5ID0gYiAtIGJ2XG4gICAgcVsxXSA9IHhcbiAgICBpZih5KSB7XG4gICAgICBnW2NvdW50KytdID0geVxuICAgIH1cbiAgfVxuICBpZihxWzFdKSB7XG4gICAgZ1tjb3VudCsrXSA9IHFbMV1cbiAgfVxuICBpZihjb3VudCA9PT0gMCkge1xuICAgIGdbY291bnQrK10gPSAwLjBcbiAgfVxuICBnLmxlbmd0aCA9IGNvdW50XG4gIHJldHVybiBnXG59IiwiXCJ1c2Ugc3RyaWN0XCJcblxubW9kdWxlLmV4cG9ydHMgPSByb2J1c3RTdWJ0cmFjdFxuXG4vL0Vhc3kgY2FzZTogQWRkIHR3byBzY2FsYXJzXG5mdW5jdGlvbiBzY2FsYXJTY2FsYXIoYSwgYikge1xuICB2YXIgeCA9IGEgKyBiXG4gIHZhciBidiA9IHggLSBhXG4gIHZhciBhdiA9IHggLSBidlxuICB2YXIgYnIgPSBiIC0gYnZcbiAgdmFyIGFyID0gYSAtIGF2XG4gIHZhciB5ID0gYXIgKyBiclxuICBpZih5KSB7XG4gICAgcmV0dXJuIFt5LCB4XVxuICB9XG4gIHJldHVybiBbeF1cbn1cblxuZnVuY3Rpb24gcm9idXN0U3VidHJhY3QoZSwgZikge1xuICB2YXIgbmUgPSBlLmxlbmd0aHwwXG4gIHZhciBuZiA9IGYubGVuZ3RofDBcbiAgaWYobmUgPT09IDEgJiYgbmYgPT09IDEpIHtcbiAgICByZXR1cm4gc2NhbGFyU2NhbGFyKGVbMF0sIC1mWzBdKVxuICB9XG4gIHZhciBuID0gbmUgKyBuZlxuICB2YXIgZyA9IG5ldyBBcnJheShuKVxuICB2YXIgY291bnQgPSAwXG4gIHZhciBlcHRyID0gMFxuICB2YXIgZnB0ciA9IDBcbiAgdmFyIGFicyA9IE1hdGguYWJzXG4gIHZhciBlaSA9IGVbZXB0cl1cbiAgdmFyIGVhID0gYWJzKGVpKVxuICB2YXIgZmkgPSAtZltmcHRyXVxuICB2YXIgZmEgPSBhYnMoZmkpXG4gIHZhciBhLCBiXG4gIGlmKGVhIDwgZmEpIHtcbiAgICBiID0gZWlcbiAgICBlcHRyICs9IDFcbiAgICBpZihlcHRyIDwgbmUpIHtcbiAgICAgIGVpID0gZVtlcHRyXVxuICAgICAgZWEgPSBhYnMoZWkpXG4gICAgfVxuICB9IGVsc2Uge1xuICAgIGIgPSBmaVxuICAgIGZwdHIgKz0gMVxuICAgIGlmKGZwdHIgPCBuZikge1xuICAgICAgZmkgPSAtZltmcHRyXVxuICAgICAgZmEgPSBhYnMoZmkpXG4gICAgfVxuICB9XG4gIGlmKChlcHRyIDwgbmUgJiYgZWEgPCBmYSkgfHwgKGZwdHIgPj0gbmYpKSB7XG4gICAgYSA9IGVpXG4gICAgZXB0ciArPSAxXG4gICAgaWYoZXB0ciA8IG5lKSB7XG4gICAgICBlaSA9IGVbZXB0cl1cbiAgICAgIGVhID0gYWJzKGVpKVxuICAgIH1cbiAgfSBlbHNlIHtcbiAgICBhID0gZmlcbiAgICBmcHRyICs9IDFcbiAgICBpZihmcHRyIDwgbmYpIHtcbiAgICAgIGZpID0gLWZbZnB0cl1cbiAgICAgIGZhID0gYWJzKGZpKVxuICAgIH1cbiAgfVxuICB2YXIgeCA9IGEgKyBiXG4gIHZhciBidiA9IHggLSBhXG4gIHZhciB5ID0gYiAtIGJ2XG4gIHZhciBxMCA9IHlcbiAgdmFyIHExID0geFxuICB2YXIgX3gsIF9idiwgX2F2LCBfYnIsIF9hclxuICB3aGlsZShlcHRyIDwgbmUgJiYgZnB0ciA8IG5mKSB7XG4gICAgaWYoZWEgPCBmYSkge1xuICAgICAgYSA9IGVpXG4gICAgICBlcHRyICs9IDFcbiAgICAgIGlmKGVwdHIgPCBuZSkge1xuICAgICAgICBlaSA9IGVbZXB0cl1cbiAgICAgICAgZWEgPSBhYnMoZWkpXG4gICAgICB9XG4gICAgfSBlbHNlIHtcbiAgICAgIGEgPSBmaVxuICAgICAgZnB0ciArPSAxXG4gICAgICBpZihmcHRyIDwgbmYpIHtcbiAgICAgICAgZmkgPSAtZltmcHRyXVxuICAgICAgICBmYSA9IGFicyhmaSlcbiAgICAgIH1cbiAgICB9XG4gICAgYiA9IHEwXG4gICAgeCA9IGEgKyBiXG4gICAgYnYgPSB4IC0gYVxuICAgIHkgPSBiIC0gYnZcbiAgICBpZih5KSB7XG4gICAgICBnW2NvdW50KytdID0geVxuICAgIH1cbiAgICBfeCA9IHExICsgeFxuICAgIF9idiA9IF94IC0gcTFcbiAgICBfYXYgPSBfeCAtIF9idlxuICAgIF9iciA9IHggLSBfYnZcbiAgICBfYXIgPSBxMSAtIF9hdlxuICAgIHEwID0gX2FyICsgX2JyXG4gICAgcTEgPSBfeFxuICB9XG4gIHdoaWxlKGVwdHIgPCBuZSkge1xuICAgIGEgPSBlaVxuICAgIGIgPSBxMFxuICAgIHggPSBhICsgYlxuICAgIGJ2ID0geCAtIGFcbiAgICB5ID0gYiAtIGJ2XG4gICAgaWYoeSkge1xuICAgICAgZ1tjb3VudCsrXSA9IHlcbiAgICB9XG4gICAgX3ggPSBxMSArIHhcbiAgICBfYnYgPSBfeCAtIHExXG4gICAgX2F2ID0gX3ggLSBfYnZcbiAgICBfYnIgPSB4IC0gX2J2XG4gICAgX2FyID0gcTEgLSBfYXZcbiAgICBxMCA9IF9hciArIF9iclxuICAgIHExID0gX3hcbiAgICBlcHRyICs9IDFcbiAgICBpZihlcHRyIDwgbmUpIHtcbiAgICAgIGVpID0gZVtlcHRyXVxuICAgIH1cbiAgfVxuICB3aGlsZShmcHRyIDwgbmYpIHtcbiAgICBhID0gZmlcbiAgICBiID0gcTBcbiAgICB4ID0gYSArIGJcbiAgICBidiA9IHggLSBhXG4gICAgeSA9IGIgLSBidlxuICAgIGlmKHkpIHtcbiAgICAgIGdbY291bnQrK10gPSB5XG4gICAgfSBcbiAgICBfeCA9IHExICsgeFxuICAgIF9idiA9IF94IC0gcTFcbiAgICBfYXYgPSBfeCAtIF9idlxuICAgIF9iciA9IHggLSBfYnZcbiAgICBfYXIgPSBxMSAtIF9hdlxuICAgIHEwID0gX2FyICsgX2JyXG4gICAgcTEgPSBfeFxuICAgIGZwdHIgKz0gMVxuICAgIGlmKGZwdHIgPCBuZikge1xuICAgICAgZmkgPSAtZltmcHRyXVxuICAgIH1cbiAgfVxuICBpZihxMCkge1xuICAgIGdbY291bnQrK10gPSBxMFxuICB9XG4gIGlmKHExKSB7XG4gICAgZ1tjb3VudCsrXSA9IHExXG4gIH1cbiAgaWYoIWNvdW50KSB7XG4gICAgZ1tjb3VudCsrXSA9IDAuMCAgXG4gIH1cbiAgZy5sZW5ndGggPSBjb3VudFxuICByZXR1cm4gZ1xufSIsIlwidXNlIHN0cmljdFwiXG5cbm1vZHVsZS5leHBvcnRzID0gbGluZWFyRXhwYW5zaW9uU3VtXG5cbi8vRWFzeSBjYXNlOiBBZGQgdHdvIHNjYWxhcnNcbmZ1bmN0aW9uIHNjYWxhclNjYWxhcihhLCBiKSB7XG4gIHZhciB4ID0gYSArIGJcbiAgdmFyIGJ2ID0geCAtIGFcbiAgdmFyIGF2ID0geCAtIGJ2XG4gIHZhciBiciA9IGIgLSBidlxuICB2YXIgYXIgPSBhIC0gYXZcbiAgdmFyIHkgPSBhciArIGJyXG4gIGlmKHkpIHtcbiAgICByZXR1cm4gW3ksIHhdXG4gIH1cbiAgcmV0dXJuIFt4XVxufVxuXG5mdW5jdGlvbiBsaW5lYXJFeHBhbnNpb25TdW0oZSwgZikge1xuICB2YXIgbmUgPSBlLmxlbmd0aHwwXG4gIHZhciBuZiA9IGYubGVuZ3RofDBcbiAgaWYobmUgPT09IDEgJiYgbmYgPT09IDEpIHtcbiAgICByZXR1cm4gc2NhbGFyU2NhbGFyKGVbMF0sIGZbMF0pXG4gIH1cbiAgdmFyIG4gPSBuZSArIG5mXG4gIHZhciBnID0gbmV3IEFycmF5KG4pXG4gIHZhciBjb3VudCA9IDBcbiAgdmFyIGVwdHIgPSAwXG4gIHZhciBmcHRyID0gMFxuICB2YXIgYWJzID0gTWF0aC5hYnNcbiAgdmFyIGVpID0gZVtlcHRyXVxuICB2YXIgZWEgPSBhYnMoZWkpXG4gIHZhciBmaSA9IGZbZnB0cl1cbiAgdmFyIGZhID0gYWJzKGZpKVxuICB2YXIgYSwgYlxuICBpZihlYSA8IGZhKSB7XG4gICAgYiA9IGVpXG4gICAgZXB0ciArPSAxXG4gICAgaWYoZXB0ciA8IG5lKSB7XG4gICAgICBlaSA9IGVbZXB0cl1cbiAgICAgIGVhID0gYWJzKGVpKVxuICAgIH1cbiAgfSBlbHNlIHtcbiAgICBiID0gZmlcbiAgICBmcHRyICs9IDFcbiAgICBpZihmcHRyIDwgbmYpIHtcbiAgICAgIGZpID0gZltmcHRyXVxuICAgICAgZmEgPSBhYnMoZmkpXG4gICAgfVxuICB9XG4gIGlmKChlcHRyIDwgbmUgJiYgZWEgPCBmYSkgfHwgKGZwdHIgPj0gbmYpKSB7XG4gICAgYSA9IGVpXG4gICAgZXB0ciArPSAxXG4gICAgaWYoZXB0ciA8IG5lKSB7XG4gICAgICBlaSA9IGVbZXB0cl1cbiAgICAgIGVhID0gYWJzKGVpKVxuICAgIH1cbiAgfSBlbHNlIHtcbiAgICBhID0gZmlcbiAgICBmcHRyICs9IDFcbiAgICBpZihmcHRyIDwgbmYpIHtcbiAgICAgIGZpID0gZltmcHRyXVxuICAgICAgZmEgPSBhYnMoZmkpXG4gICAgfVxuICB9XG4gIHZhciB4ID0gYSArIGJcbiAgdmFyIGJ2ID0geCAtIGFcbiAgdmFyIHkgPSBiIC0gYnZcbiAgdmFyIHEwID0geVxuICB2YXIgcTEgPSB4XG4gIHZhciBfeCwgX2J2LCBfYXYsIF9iciwgX2FyXG4gIHdoaWxlKGVwdHIgPCBuZSAmJiBmcHRyIDwgbmYpIHtcbiAgICBpZihlYSA8IGZhKSB7XG4gICAgICBhID0gZWlcbiAgICAgIGVwdHIgKz0gMVxuICAgICAgaWYoZXB0ciA8IG5lKSB7XG4gICAgICAgIGVpID0gZVtlcHRyXVxuICAgICAgICBlYSA9IGFicyhlaSlcbiAgICAgIH1cbiAgICB9IGVsc2Uge1xuICAgICAgYSA9IGZpXG4gICAgICBmcHRyICs9IDFcbiAgICAgIGlmKGZwdHIgPCBuZikge1xuICAgICAgICBmaSA9IGZbZnB0cl1cbiAgICAgICAgZmEgPSBhYnMoZmkpXG4gICAgICB9XG4gICAgfVxuICAgIGIgPSBxMFxuICAgIHggPSBhICsgYlxuICAgIGJ2ID0geCAtIGFcbiAgICB5ID0gYiAtIGJ2XG4gICAgaWYoeSkge1xuICAgICAgZ1tjb3VudCsrXSA9IHlcbiAgICB9XG4gICAgX3ggPSBxMSArIHhcbiAgICBfYnYgPSBfeCAtIHExXG4gICAgX2F2ID0gX3ggLSBfYnZcbiAgICBfYnIgPSB4IC0gX2J2XG4gICAgX2FyID0gcTEgLSBfYXZcbiAgICBxMCA9IF9hciArIF9iclxuICAgIHExID0gX3hcbiAgfVxuICB3aGlsZShlcHRyIDwgbmUpIHtcbiAgICBhID0gZWlcbiAgICBiID0gcTBcbiAgICB4ID0gYSArIGJcbiAgICBidiA9IHggLSBhXG4gICAgeSA9IGIgLSBidlxuICAgIGlmKHkpIHtcbiAgICAgIGdbY291bnQrK10gPSB5XG4gICAgfVxuICAgIF94ID0gcTEgKyB4XG4gICAgX2J2ID0gX3ggLSBxMVxuICAgIF9hdiA9IF94IC0gX2J2XG4gICAgX2JyID0geCAtIF9idlxuICAgIF9hciA9IHExIC0gX2F2XG4gICAgcTAgPSBfYXIgKyBfYnJcbiAgICBxMSA9IF94XG4gICAgZXB0ciArPSAxXG4gICAgaWYoZXB0ciA8IG5lKSB7XG4gICAgICBlaSA9IGVbZXB0cl1cbiAgICB9XG4gIH1cbiAgd2hpbGUoZnB0ciA8IG5mKSB7XG4gICAgYSA9IGZpXG4gICAgYiA9IHEwXG4gICAgeCA9IGEgKyBiXG4gICAgYnYgPSB4IC0gYVxuICAgIHkgPSBiIC0gYnZcbiAgICBpZih5KSB7XG4gICAgICBnW2NvdW50KytdID0geVxuICAgIH0gXG4gICAgX3ggPSBxMSArIHhcbiAgICBfYnYgPSBfeCAtIHExXG4gICAgX2F2ID0gX3ggLSBfYnZcbiAgICBfYnIgPSB4IC0gX2J2XG4gICAgX2FyID0gcTEgLSBfYXZcbiAgICBxMCA9IF9hciArIF9iclxuICAgIHExID0gX3hcbiAgICBmcHRyICs9IDFcbiAgICBpZihmcHRyIDwgbmYpIHtcbiAgICAgIGZpID0gZltmcHRyXVxuICAgIH1cbiAgfVxuICBpZihxMCkge1xuICAgIGdbY291bnQrK10gPSBxMFxuICB9XG4gIGlmKHExKSB7XG4gICAgZ1tjb3VudCsrXSA9IHExXG4gIH1cbiAgaWYoIWNvdW50KSB7XG4gICAgZ1tjb3VudCsrXSA9IDAuMCAgXG4gIH1cbiAgZy5sZW5ndGggPSBjb3VudFxuICByZXR1cm4gZ1xufSIsIlwidXNlIHN0cmljdFwiXG5cbm1vZHVsZS5leHBvcnRzID0gdHdvUHJvZHVjdFxuXG52YXIgU1BMSVRURVIgPSArKE1hdGgucG93KDIsIDI3KSArIDEuMClcblxuZnVuY3Rpb24gdHdvUHJvZHVjdChhLCBiLCByZXN1bHQpIHtcbiAgdmFyIHggPSBhICogYlxuXG4gIHZhciBjID0gU1BMSVRURVIgKiBhXG4gIHZhciBhYmlnID0gYyAtIGFcbiAgdmFyIGFoaSA9IGMgLSBhYmlnXG4gIHZhciBhbG8gPSBhIC0gYWhpXG5cbiAgdmFyIGQgPSBTUExJVFRFUiAqIGJcbiAgdmFyIGJiaWcgPSBkIC0gYlxuICB2YXIgYmhpID0gZCAtIGJiaWdcbiAgdmFyIGJsbyA9IGIgLSBiaGlcblxuICB2YXIgZXJyMSA9IHggLSAoYWhpICogYmhpKVxuICB2YXIgZXJyMiA9IGVycjEgLSAoYWxvICogYmhpKVxuICB2YXIgZXJyMyA9IGVycjIgLSAoYWhpICogYmxvKVxuXG4gIHZhciB5ID0gYWxvICogYmxvIC0gZXJyM1xuXG4gIGlmKHJlc3VsdCkge1xuICAgIHJlc3VsdFswXSA9IHlcbiAgICByZXN1bHRbMV0gPSB4XG4gICAgcmV0dXJuIHJlc3VsdFxuICB9XG5cbiAgcmV0dXJuIFsgeSwgeCBdXG59IiwiXCJ1c2Ugc3RyaWN0XCJcblxudmFyIHR3b1Byb2R1Y3QgPSByZXF1aXJlKFwidHdvLXByb2R1Y3RcIilcbnZhciByb2J1c3RTdW0gPSByZXF1aXJlKFwicm9idXN0LXN1bVwiKVxudmFyIHJvYnVzdFNjYWxlID0gcmVxdWlyZShcInJvYnVzdC1zY2FsZVwiKVxudmFyIHJvYnVzdFN1YnRyYWN0ID0gcmVxdWlyZShcInJvYnVzdC1zdWJ0cmFjdFwiKVxuXG52YXIgTlVNX0VYUEFORCA9IDVcblxudmFyIEVQU0lMT04gICAgID0gMS4xMTAyMjMwMjQ2MjUxNTY1ZS0xNlxudmFyIEVSUkJPVU5EMyAgID0gKDMuMCArIDE2LjAgKiBFUFNJTE9OKSAqIEVQU0lMT05cbnZhciBFUlJCT1VORDQgICA9ICg3LjAgKyA1Ni4wICogRVBTSUxPTikgKiBFUFNJTE9OXG5cbmZ1bmN0aW9uIGNvZmFjdG9yKG0sIGMpIHtcbiAgdmFyIHJlc3VsdCA9IG5ldyBBcnJheShtLmxlbmd0aC0xKVxuICBmb3IodmFyIGk9MTsgaTxtLmxlbmd0aDsgKytpKSB7XG4gICAgdmFyIHIgPSByZXN1bHRbaS0xXSA9IG5ldyBBcnJheShtLmxlbmd0aC0xKVxuICAgIGZvcih2YXIgaj0wLGs9MDsgajxtLmxlbmd0aDsgKytqKSB7XG4gICAgICBpZihqID09PSBjKSB7XG4gICAgICAgIGNvbnRpbnVlXG4gICAgICB9XG4gICAgICByW2srK10gPSBtW2ldW2pdXG4gICAgfVxuICB9XG4gIHJldHVybiByZXN1bHRcbn1cblxuZnVuY3Rpb24gbWF0cml4KG4pIHtcbiAgdmFyIHJlc3VsdCA9IG5ldyBBcnJheShuKVxuICBmb3IodmFyIGk9MDsgaTxuOyArK2kpIHtcbiAgICByZXN1bHRbaV0gPSBuZXcgQXJyYXkobilcbiAgICBmb3IodmFyIGo9MDsgajxuOyArK2opIHtcbiAgICAgIHJlc3VsdFtpXVtqXSA9IFtcIm1cIiwgaiwgXCJbXCIsIChuLWktMSksIFwiXVwiXS5qb2luKFwiXCIpXG4gICAgfVxuICB9XG4gIHJldHVybiByZXN1bHRcbn1cblxuZnVuY3Rpb24gc2lnbihuKSB7XG4gIGlmKG4gJiAxKSB7XG4gICAgcmV0dXJuIFwiLVwiXG4gIH1cbiAgcmV0dXJuIFwiXCJcbn1cblxuZnVuY3Rpb24gZ2VuZXJhdGVTdW0oZXhwcikge1xuICBpZihleHByLmxlbmd0aCA9PT0gMSkge1xuICAgIHJldHVybiBleHByWzBdXG4gIH0gZWxzZSBpZihleHByLmxlbmd0aCA9PT0gMikge1xuICAgIHJldHVybiBbXCJzdW0oXCIsIGV4cHJbMF0sIFwiLFwiLCBleHByWzFdLCBcIilcIl0uam9pbihcIlwiKVxuICB9IGVsc2Uge1xuICAgIHZhciBtID0gZXhwci5sZW5ndGg+PjFcbiAgICByZXR1cm4gW1wic3VtKFwiLCBnZW5lcmF0ZVN1bShleHByLnNsaWNlKDAsIG0pKSwgXCIsXCIsIGdlbmVyYXRlU3VtKGV4cHIuc2xpY2UobSkpLCBcIilcIl0uam9pbihcIlwiKVxuICB9XG59XG5cbmZ1bmN0aW9uIGRldGVybWluYW50KG0pIHtcbiAgaWYobS5sZW5ndGggPT09IDIpIHtcbiAgICByZXR1cm4gW1tcInN1bShwcm9kKFwiLCBtWzBdWzBdLCBcIixcIiwgbVsxXVsxXSwgXCIpLHByb2QoLVwiLCBtWzBdWzFdLCBcIixcIiwgbVsxXVswXSwgXCIpKVwiXS5qb2luKFwiXCIpXVxuICB9IGVsc2Uge1xuICAgIHZhciBleHByID0gW11cbiAgICBmb3IodmFyIGk9MDsgaTxtLmxlbmd0aDsgKytpKSB7XG4gICAgICBleHByLnB1c2goW1wic2NhbGUoXCIsIGdlbmVyYXRlU3VtKGRldGVybWluYW50KGNvZmFjdG9yKG0sIGkpKSksIFwiLFwiLCBzaWduKGkpLCBtWzBdW2ldLCBcIilcIl0uam9pbihcIlwiKSlcbiAgICB9XG4gICAgcmV0dXJuIGV4cHJcbiAgfVxufVxuXG5mdW5jdGlvbiBvcmllbnRhdGlvbihuKSB7XG4gIHZhciBwb3MgPSBbXVxuICB2YXIgbmVnID0gW11cbiAgdmFyIG0gPSBtYXRyaXgobilcbiAgdmFyIGFyZ3MgPSBbXVxuICBmb3IodmFyIGk9MDsgaTxuOyArK2kpIHtcbiAgICBpZigoaSYxKT09PTApIHtcbiAgICAgIHBvcy5wdXNoLmFwcGx5KHBvcywgZGV0ZXJtaW5hbnQoY29mYWN0b3IobSwgaSkpKVxuICAgIH0gZWxzZSB7XG4gICAgICBuZWcucHVzaC5hcHBseShuZWcsIGRldGVybWluYW50KGNvZmFjdG9yKG0sIGkpKSlcbiAgICB9XG4gICAgYXJncy5wdXNoKFwibVwiICsgaSlcbiAgfVxuICB2YXIgcG9zRXhwciA9IGdlbmVyYXRlU3VtKHBvcylcbiAgdmFyIG5lZ0V4cHIgPSBnZW5lcmF0ZVN1bShuZWcpXG4gIHZhciBmdW5jTmFtZSA9IFwib3JpZW50YXRpb25cIiArIG4gKyBcIkV4YWN0XCJcbiAgdmFyIGNvZGUgPSBbXCJmdW5jdGlvbiBcIiwgZnVuY05hbWUsIFwiKFwiLCBhcmdzLmpvaW4oKSwgXCIpe3ZhciBwPVwiLCBwb3NFeHByLCBcIixuPVwiLCBuZWdFeHByLCBcIixkPXN1YihwLG4pO1xcXG5yZXR1cm4gZFtkLmxlbmd0aC0xXTt9O3JldHVybiBcIiwgZnVuY05hbWVdLmpvaW4oXCJcIilcbiAgdmFyIHByb2MgPSBuZXcgRnVuY3Rpb24oXCJzdW1cIiwgXCJwcm9kXCIsIFwic2NhbGVcIiwgXCJzdWJcIiwgY29kZSlcbiAgcmV0dXJuIHByb2Mocm9idXN0U3VtLCB0d29Qcm9kdWN0LCByb2J1c3RTY2FsZSwgcm9idXN0U3VidHJhY3QpXG59XG5cbnZhciBvcmllbnRhdGlvbjNFeGFjdCA9IG9yaWVudGF0aW9uKDMpXG52YXIgb3JpZW50YXRpb240RXhhY3QgPSBvcmllbnRhdGlvbig0KVxuXG52YXIgQ0FDSEVEID0gW1xuICBmdW5jdGlvbiBvcmllbnRhdGlvbjAoKSB7IHJldHVybiAwIH0sXG4gIGZ1bmN0aW9uIG9yaWVudGF0aW9uMSgpIHsgcmV0dXJuIDAgfSxcbiAgZnVuY3Rpb24gb3JpZW50YXRpb24yKGEsIGIpIHsgXG4gICAgcmV0dXJuIGJbMF0gLSBhWzBdXG4gIH0sXG4gIGZ1bmN0aW9uIG9yaWVudGF0aW9uMyhhLCBiLCBjKSB7XG4gICAgdmFyIGwgPSAoYVsxXSAtIGNbMV0pICogKGJbMF0gLSBjWzBdKVxuICAgIHZhciByID0gKGFbMF0gLSBjWzBdKSAqIChiWzFdIC0gY1sxXSlcbiAgICB2YXIgZGV0ID0gbCAtIHJcbiAgICB2YXIgc1xuICAgIGlmKGwgPiAwKSB7XG4gICAgICBpZihyIDw9IDApIHtcbiAgICAgICAgcmV0dXJuIGRldFxuICAgICAgfSBlbHNlIHtcbiAgICAgICAgcyA9IGwgKyByXG4gICAgICB9XG4gICAgfSBlbHNlIGlmKGwgPCAwKSB7XG4gICAgICBpZihyID49IDApIHtcbiAgICAgICAgcmV0dXJuIGRldFxuICAgICAgfSBlbHNlIHtcbiAgICAgICAgcyA9IC0obCArIHIpXG4gICAgICB9XG4gICAgfSBlbHNlIHtcbiAgICAgIHJldHVybiBkZXRcbiAgICB9XG4gICAgdmFyIHRvbCA9IEVSUkJPVU5EMyAqIHNcbiAgICBpZihkZXQgPj0gdG9sIHx8IGRldCA8PSAtdG9sKSB7XG4gICAgICByZXR1cm4gZGV0XG4gICAgfVxuICAgIHJldHVybiBvcmllbnRhdGlvbjNFeGFjdChhLCBiLCBjKVxuICB9LFxuICBmdW5jdGlvbiBvcmllbnRhdGlvbjQoYSxiLGMsZCkge1xuICAgIHZhciBhZHggPSBhWzBdIC0gZFswXVxuICAgIHZhciBiZHggPSBiWzBdIC0gZFswXVxuICAgIHZhciBjZHggPSBjWzBdIC0gZFswXVxuICAgIHZhciBhZHkgPSBhWzFdIC0gZFsxXVxuICAgIHZhciBiZHkgPSBiWzFdIC0gZFsxXVxuICAgIHZhciBjZHkgPSBjWzFdIC0gZFsxXVxuICAgIHZhciBhZHogPSBhWzJdIC0gZFsyXVxuICAgIHZhciBiZHogPSBiWzJdIC0gZFsyXVxuICAgIHZhciBjZHogPSBjWzJdIC0gZFsyXVxuICAgIHZhciBiZHhjZHkgPSBiZHggKiBjZHlcbiAgICB2YXIgY2R4YmR5ID0gY2R4ICogYmR5XG4gICAgdmFyIGNkeGFkeSA9IGNkeCAqIGFkeVxuICAgIHZhciBhZHhjZHkgPSBhZHggKiBjZHlcbiAgICB2YXIgYWR4YmR5ID0gYWR4ICogYmR5XG4gICAgdmFyIGJkeGFkeSA9IGJkeCAqIGFkeVxuICAgIHZhciBkZXQgPSBhZHogKiAoYmR4Y2R5IC0gY2R4YmR5KSBcbiAgICAgICAgICAgICsgYmR6ICogKGNkeGFkeSAtIGFkeGNkeSlcbiAgICAgICAgICAgICsgY2R6ICogKGFkeGJkeSAtIGJkeGFkeSlcbiAgICB2YXIgcGVybWFuZW50ID0gKE1hdGguYWJzKGJkeGNkeSkgKyBNYXRoLmFicyhjZHhiZHkpKSAqIE1hdGguYWJzKGFkeilcbiAgICAgICAgICAgICAgICAgICsgKE1hdGguYWJzKGNkeGFkeSkgKyBNYXRoLmFicyhhZHhjZHkpKSAqIE1hdGguYWJzKGJkeilcbiAgICAgICAgICAgICAgICAgICsgKE1hdGguYWJzKGFkeGJkeSkgKyBNYXRoLmFicyhiZHhhZHkpKSAqIE1hdGguYWJzKGNkeilcbiAgICB2YXIgdG9sID0gRVJSQk9VTkQ0ICogcGVybWFuZW50XG4gICAgaWYgKChkZXQgPiB0b2wpIHx8ICgtZGV0ID4gdG9sKSkge1xuICAgICAgcmV0dXJuIGRldFxuICAgIH1cbiAgICByZXR1cm4gb3JpZW50YXRpb240RXhhY3QoYSxiLGMsZClcbiAgfVxuXVxuXG5mdW5jdGlvbiBzbG93T3JpZW50KGFyZ3MpIHtcbiAgdmFyIHByb2MgPSBDQUNIRURbYXJncy5sZW5ndGhdXG4gIGlmKCFwcm9jKSB7XG4gICAgcHJvYyA9IENBQ0hFRFthcmdzLmxlbmd0aF0gPSBvcmllbnRhdGlvbihhcmdzLmxlbmd0aClcbiAgfVxuICByZXR1cm4gcHJvYy5hcHBseSh1bmRlZmluZWQsIGFyZ3MpXG59XG5cbmZ1bmN0aW9uIGdlbmVyYXRlT3JpZW50YXRpb25Qcm9jKCkge1xuICB3aGlsZShDQUNIRUQubGVuZ3RoIDw9IE5VTV9FWFBBTkQpIHtcbiAgICBDQUNIRUQucHVzaChvcmllbnRhdGlvbihDQUNIRUQubGVuZ3RoKSlcbiAgfVxuICB2YXIgYXJncyA9IFtdXG4gIHZhciBwcm9jQXJncyA9IFtcInNsb3dcIl1cbiAgZm9yKHZhciBpPTA7IGk8PU5VTV9FWFBBTkQ7ICsraSkge1xuICAgIGFyZ3MucHVzaChcImFcIiArIGkpXG4gICAgcHJvY0FyZ3MucHVzaChcIm9cIiArIGkpXG4gIH1cbiAgdmFyIGNvZGUgPSBbXG4gICAgXCJmdW5jdGlvbiBnZXRPcmllbnRhdGlvbihcIiwgYXJncy5qb2luKCksIFwiKXtzd2l0Y2goYXJndW1lbnRzLmxlbmd0aCl7Y2FzZSAwOmNhc2UgMTpyZXR1cm4gMDtcIlxuICBdXG4gIGZvcih2YXIgaT0yOyBpPD1OVU1fRVhQQU5EOyArK2kpIHtcbiAgICBjb2RlLnB1c2goXCJjYXNlIFwiLCBpLCBcIjpyZXR1cm4gb1wiLCBpLCBcIihcIiwgYXJncy5zbGljZSgwLCBpKS5qb2luKCksIFwiKTtcIilcbiAgfVxuICBjb2RlLnB1c2goXCJ9dmFyIHM9bmV3IEFycmF5KGFyZ3VtZW50cy5sZW5ndGgpO2Zvcih2YXIgaT0wO2k8YXJndW1lbnRzLmxlbmd0aDsrK2kpe3NbaV09YXJndW1lbnRzW2ldfTtyZXR1cm4gc2xvdyhzKTt9cmV0dXJuIGdldE9yaWVudGF0aW9uXCIpXG4gIHByb2NBcmdzLnB1c2goY29kZS5qb2luKFwiXCIpKVxuXG4gIHZhciBwcm9jID0gRnVuY3Rpb24uYXBwbHkodW5kZWZpbmVkLCBwcm9jQXJncylcbiAgbW9kdWxlLmV4cG9ydHMgPSBwcm9jLmFwcGx5KHVuZGVmaW5lZCwgW3Nsb3dPcmllbnRdLmNvbmNhdChDQUNIRUQpKVxuICBmb3IodmFyIGk9MDsgaTw9TlVNX0VYUEFORDsgKytpKSB7XG4gICAgbW9kdWxlLmV4cG9ydHNbaV0gPSBDQUNIRURbaV1cbiAgfVxufVxuXG5nZW5lcmF0ZU9yaWVudGF0aW9uUHJvYygpIiwiLyoqXG4gKiBCaXQgdHdpZGRsaW5nIGhhY2tzIGZvciBKYXZhU2NyaXB0LlxuICpcbiAqIEF1dGhvcjogTWlrb2xhIEx5c2Vua29cbiAqXG4gKiBQb3J0ZWQgZnJvbSBTdGFuZm9yZCBiaXQgdHdpZGRsaW5nIGhhY2sgbGlicmFyeTpcbiAqICAgIGh0dHA6Ly9ncmFwaGljcy5zdGFuZm9yZC5lZHUvfnNlYW5kZXIvYml0aGFja3MuaHRtbFxuICovXG5cblwidXNlIHN0cmljdFwiOyBcInVzZSByZXN0cmljdFwiO1xuXG4vL051bWJlciBvZiBiaXRzIGluIGFuIGludGVnZXJcbnZhciBJTlRfQklUUyA9IDMyO1xuXG4vL0NvbnN0YW50c1xuZXhwb3J0cy5JTlRfQklUUyAgPSBJTlRfQklUUztcbmV4cG9ydHMuSU5UX01BWCAgID0gIDB4N2ZmZmZmZmY7XG5leHBvcnRzLklOVF9NSU4gICA9IC0xPDwoSU5UX0JJVFMtMSk7XG5cbi8vUmV0dXJucyAtMSwgMCwgKzEgZGVwZW5kaW5nIG9uIHNpZ24gb2YgeFxuZXhwb3J0cy5zaWduID0gZnVuY3Rpb24odikge1xuICByZXR1cm4gKHYgPiAwKSAtICh2IDwgMCk7XG59XG5cbi8vQ29tcHV0ZXMgYWJzb2x1dGUgdmFsdWUgb2YgaW50ZWdlclxuZXhwb3J0cy5hYnMgPSBmdW5jdGlvbih2KSB7XG4gIHZhciBtYXNrID0gdiA+PiAoSU5UX0JJVFMtMSk7XG4gIHJldHVybiAodiBeIG1hc2spIC0gbWFzaztcbn1cblxuLy9Db21wdXRlcyBtaW5pbXVtIG9mIGludGVnZXJzIHggYW5kIHlcbmV4cG9ydHMubWluID0gZnVuY3Rpb24oeCwgeSkge1xuICByZXR1cm4geSBeICgoeCBeIHkpICYgLSh4IDwgeSkpO1xufVxuXG4vL0NvbXB1dGVzIG1heGltdW0gb2YgaW50ZWdlcnMgeCBhbmQgeVxuZXhwb3J0cy5tYXggPSBmdW5jdGlvbih4LCB5KSB7XG4gIHJldHVybiB4IF4gKCh4IF4geSkgJiAtKHggPCB5KSk7XG59XG5cbi8vQ2hlY2tzIGlmIGEgbnVtYmVyIGlzIGEgcG93ZXIgb2YgdHdvXG5leHBvcnRzLmlzUG93MiA9IGZ1bmN0aW9uKHYpIHtcbiAgcmV0dXJuICEodiAmICh2LTEpKSAmJiAoISF2KTtcbn1cblxuLy9Db21wdXRlcyBsb2cgYmFzZSAyIG9mIHZcbmV4cG9ydHMubG9nMiA9IGZ1bmN0aW9uKHYpIHtcbiAgdmFyIHIsIHNoaWZ0O1xuICByID0gICAgICh2ID4gMHhGRkZGKSA8PCA0OyB2ID4+Pj0gcjtcbiAgc2hpZnQgPSAodiA+IDB4RkYgICkgPDwgMzsgdiA+Pj49IHNoaWZ0OyByIHw9IHNoaWZ0O1xuICBzaGlmdCA9ICh2ID4gMHhGICAgKSA8PCAyOyB2ID4+Pj0gc2hpZnQ7IHIgfD0gc2hpZnQ7XG4gIHNoaWZ0ID0gKHYgPiAweDMgICApIDw8IDE7IHYgPj4+PSBzaGlmdDsgciB8PSBzaGlmdDtcbiAgcmV0dXJuIHIgfCAodiA+PiAxKTtcbn1cblxuLy9Db21wdXRlcyBsb2cgYmFzZSAxMCBvZiB2XG5leHBvcnRzLmxvZzEwID0gZnVuY3Rpb24odikge1xuICByZXR1cm4gICh2ID49IDEwMDAwMDAwMDApID8gOSA6ICh2ID49IDEwMDAwMDAwMCkgPyA4IDogKHYgPj0gMTAwMDAwMDApID8gNyA6XG4gICAgICAgICAgKHYgPj0gMTAwMDAwMCkgPyA2IDogKHYgPj0gMTAwMDAwKSA/IDUgOiAodiA+PSAxMDAwMCkgPyA0IDpcbiAgICAgICAgICAodiA+PSAxMDAwKSA/IDMgOiAodiA+PSAxMDApID8gMiA6ICh2ID49IDEwKSA/IDEgOiAwO1xufVxuXG4vL0NvdW50cyBudW1iZXIgb2YgYml0c1xuZXhwb3J0cy5wb3BDb3VudCA9IGZ1bmN0aW9uKHYpIHtcbiAgdiA9IHYgLSAoKHYgPj4+IDEpICYgMHg1NTU1NTU1NSk7XG4gIHYgPSAodiAmIDB4MzMzMzMzMzMpICsgKCh2ID4+PiAyKSAmIDB4MzMzMzMzMzMpO1xuICByZXR1cm4gKCh2ICsgKHYgPj4+IDQpICYgMHhGMEYwRjBGKSAqIDB4MTAxMDEwMSkgPj4+IDI0O1xufVxuXG4vL0NvdW50cyBudW1iZXIgb2YgdHJhaWxpbmcgemVyb3NcbmZ1bmN0aW9uIGNvdW50VHJhaWxpbmdaZXJvcyh2KSB7XG4gIHZhciBjID0gMzI7XG4gIHYgJj0gLXY7XG4gIGlmICh2KSBjLS07XG4gIGlmICh2ICYgMHgwMDAwRkZGRikgYyAtPSAxNjtcbiAgaWYgKHYgJiAweDAwRkYwMEZGKSBjIC09IDg7XG4gIGlmICh2ICYgMHgwRjBGMEYwRikgYyAtPSA0O1xuICBpZiAodiAmIDB4MzMzMzMzMzMpIGMgLT0gMjtcbiAgaWYgKHYgJiAweDU1NTU1NTU1KSBjIC09IDE7XG4gIHJldHVybiBjO1xufVxuZXhwb3J0cy5jb3VudFRyYWlsaW5nWmVyb3MgPSBjb3VudFRyYWlsaW5nWmVyb3M7XG5cbi8vUm91bmRzIHRvIG5leHQgcG93ZXIgb2YgMlxuZXhwb3J0cy5uZXh0UG93MiA9IGZ1bmN0aW9uKHYpIHtcbiAgdiArPSB2ID09PSAwO1xuICAtLXY7XG4gIHYgfD0gdiA+Pj4gMTtcbiAgdiB8PSB2ID4+PiAyO1xuICB2IHw9IHYgPj4+IDQ7XG4gIHYgfD0gdiA+Pj4gODtcbiAgdiB8PSB2ID4+PiAxNjtcbiAgcmV0dXJuIHYgKyAxO1xufVxuXG4vL1JvdW5kcyBkb3duIHRvIHByZXZpb3VzIHBvd2VyIG9mIDJcbmV4cG9ydHMucHJldlBvdzIgPSBmdW5jdGlvbih2KSB7XG4gIHYgfD0gdiA+Pj4gMTtcbiAgdiB8PSB2ID4+PiAyO1xuICB2IHw9IHYgPj4+IDQ7XG4gIHYgfD0gdiA+Pj4gODtcbiAgdiB8PSB2ID4+PiAxNjtcbiAgcmV0dXJuIHYgLSAodj4+PjEpO1xufVxuXG4vL0NvbXB1dGVzIHBhcml0eSBvZiB3b3JkXG5leHBvcnRzLnBhcml0eSA9IGZ1bmN0aW9uKHYpIHtcbiAgdiBePSB2ID4+PiAxNjtcbiAgdiBePSB2ID4+PiA4O1xuICB2IF49IHYgPj4+IDQ7XG4gIHYgJj0gMHhmO1xuICByZXR1cm4gKDB4Njk5NiA+Pj4gdikgJiAxO1xufVxuXG52YXIgUkVWRVJTRV9UQUJMRSA9IG5ldyBBcnJheSgyNTYpO1xuXG4oZnVuY3Rpb24odGFiKSB7XG4gIGZvcih2YXIgaT0wOyBpPDI1NjsgKytpKSB7XG4gICAgdmFyIHYgPSBpLCByID0gaSwgcyA9IDc7XG4gICAgZm9yICh2ID4+Pj0gMTsgdjsgdiA+Pj49IDEpIHtcbiAgICAgIHIgPDw9IDE7XG4gICAgICByIHw9IHYgJiAxO1xuICAgICAgLS1zO1xuICAgIH1cbiAgICB0YWJbaV0gPSAociA8PCBzKSAmIDB4ZmY7XG4gIH1cbn0pKFJFVkVSU0VfVEFCTEUpO1xuXG4vL1JldmVyc2UgYml0cyBpbiBhIDMyIGJpdCB3b3JkXG5leHBvcnRzLnJldmVyc2UgPSBmdW5jdGlvbih2KSB7XG4gIHJldHVybiAgKFJFVkVSU0VfVEFCTEVbIHYgICAgICAgICAmIDB4ZmZdIDw8IDI0KSB8XG4gICAgICAgICAgKFJFVkVSU0VfVEFCTEVbKHYgPj4+IDgpICAmIDB4ZmZdIDw8IDE2KSB8XG4gICAgICAgICAgKFJFVkVSU0VfVEFCTEVbKHYgPj4+IDE2KSAmIDB4ZmZdIDw8IDgpICB8XG4gICAgICAgICAgIFJFVkVSU0VfVEFCTEVbKHYgPj4+IDI0KSAmIDB4ZmZdO1xufVxuXG4vL0ludGVybGVhdmUgYml0cyBvZiAyIGNvb3JkaW5hdGVzIHdpdGggMTYgYml0cy4gIFVzZWZ1bCBmb3IgZmFzdCBxdWFkdHJlZSBjb2Rlc1xuZXhwb3J0cy5pbnRlcmxlYXZlMiA9IGZ1bmN0aW9uKHgsIHkpIHtcbiAgeCAmPSAweEZGRkY7XG4gIHggPSAoeCB8ICh4IDw8IDgpKSAmIDB4MDBGRjAwRkY7XG4gIHggPSAoeCB8ICh4IDw8IDQpKSAmIDB4MEYwRjBGMEY7XG4gIHggPSAoeCB8ICh4IDw8IDIpKSAmIDB4MzMzMzMzMzM7XG4gIHggPSAoeCB8ICh4IDw8IDEpKSAmIDB4NTU1NTU1NTU7XG5cbiAgeSAmPSAweEZGRkY7XG4gIHkgPSAoeSB8ICh5IDw8IDgpKSAmIDB4MDBGRjAwRkY7XG4gIHkgPSAoeSB8ICh5IDw8IDQpKSAmIDB4MEYwRjBGMEY7XG4gIHkgPSAoeSB8ICh5IDw8IDIpKSAmIDB4MzMzMzMzMzM7XG4gIHkgPSAoeSB8ICh5IDw8IDEpKSAmIDB4NTU1NTU1NTU7XG5cbiAgcmV0dXJuIHggfCAoeSA8PCAxKTtcbn1cblxuLy9FeHRyYWN0cyB0aGUgbnRoIGludGVybGVhdmVkIGNvbXBvbmVudFxuZXhwb3J0cy5kZWludGVybGVhdmUyID0gZnVuY3Rpb24odiwgbikge1xuICB2ID0gKHYgPj4+IG4pICYgMHg1NTU1NTU1NTtcbiAgdiA9ICh2IHwgKHYgPj4+IDEpKSAgJiAweDMzMzMzMzMzO1xuICB2ID0gKHYgfCAodiA+Pj4gMikpICAmIDB4MEYwRjBGMEY7XG4gIHYgPSAodiB8ICh2ID4+PiA0KSkgICYgMHgwMEZGMDBGRjtcbiAgdiA9ICh2IHwgKHYgPj4+IDE2KSkgJiAweDAwMEZGRkY7XG4gIHJldHVybiAodiA8PCAxNikgPj4gMTY7XG59XG5cblxuLy9JbnRlcmxlYXZlIGJpdHMgb2YgMyBjb29yZGluYXRlcywgZWFjaCB3aXRoIDEwIGJpdHMuICBVc2VmdWwgZm9yIGZhc3Qgb2N0cmVlIGNvZGVzXG5leHBvcnRzLmludGVybGVhdmUzID0gZnVuY3Rpb24oeCwgeSwgeikge1xuICB4ICY9IDB4M0ZGO1xuICB4ICA9ICh4IHwgKHg8PDE2KSkgJiA0Mjc4MTkwMzM1O1xuICB4ICA9ICh4IHwgKHg8PDgpKSAgJiAyNTE3MTk2OTU7XG4gIHggID0gKHggfCAoeDw8NCkpICAmIDMyNzIzNTYwMzU7XG4gIHggID0gKHggfCAoeDw8MikpICAmIDEyMjcxMzM1MTM7XG5cbiAgeSAmPSAweDNGRjtcbiAgeSAgPSAoeSB8ICh5PDwxNikpICYgNDI3ODE5MDMzNTtcbiAgeSAgPSAoeSB8ICh5PDw4KSkgICYgMjUxNzE5Njk1O1xuICB5ICA9ICh5IHwgKHk8PDQpKSAgJiAzMjcyMzU2MDM1O1xuICB5ICA9ICh5IHwgKHk8PDIpKSAgJiAxMjI3MTMzNTEzO1xuICB4IHw9ICh5IDw8IDEpO1xuICBcbiAgeiAmPSAweDNGRjtcbiAgeiAgPSAoeiB8ICh6PDwxNikpICYgNDI3ODE5MDMzNTtcbiAgeiAgPSAoeiB8ICh6PDw4KSkgICYgMjUxNzE5Njk1O1xuICB6ICA9ICh6IHwgKHo8PDQpKSAgJiAzMjcyMzU2MDM1O1xuICB6ICA9ICh6IHwgKHo8PDIpKSAgJiAxMjI3MTMzNTEzO1xuICBcbiAgcmV0dXJuIHggfCAoeiA8PCAyKTtcbn1cblxuLy9FeHRyYWN0cyBudGggaW50ZXJsZWF2ZWQgY29tcG9uZW50IG9mIGEgMy10dXBsZVxuZXhwb3J0cy5kZWludGVybGVhdmUzID0gZnVuY3Rpb24odiwgbikge1xuICB2ID0gKHYgPj4+IG4pICAgICAgICYgMTIyNzEzMzUxMztcbiAgdiA9ICh2IHwgKHY+Pj4yKSkgICAmIDMyNzIzNTYwMzU7XG4gIHYgPSAodiB8ICh2Pj4+NCkpICAgJiAyNTE3MTk2OTU7XG4gIHYgPSAodiB8ICh2Pj4+OCkpICAgJiA0Mjc4MTkwMzM1O1xuICB2ID0gKHYgfCAodj4+PjE2KSkgICYgMHgzRkY7XG4gIHJldHVybiAodjw8MjIpPj4yMjtcbn1cblxuLy9Db21wdXRlcyBuZXh0IGNvbWJpbmF0aW9uIGluIGNvbGV4aWNvZ3JhcGhpYyBvcmRlciAodGhpcyBpcyBtaXN0YWtlbmx5IGNhbGxlZCBuZXh0UGVybXV0YXRpb24gb24gdGhlIGJpdCB0d2lkZGxpbmcgaGFja3MgcGFnZSlcbmV4cG9ydHMubmV4dENvbWJpbmF0aW9uID0gZnVuY3Rpb24odikge1xuICB2YXIgdCA9IHYgfCAodiAtIDEpO1xuICByZXR1cm4gKHQgKyAxKSB8ICgoKH50ICYgLX50KSAtIDEpID4+PiAoY291bnRUcmFpbGluZ1plcm9zKHYpICsgMSkpO1xufVxuXG4iLCJcInVzZSBzdHJpY3RcIjsgXCJ1c2UgcmVzdHJpY3RcIjtcblxubW9kdWxlLmV4cG9ydHMgPSBVbmlvbkZpbmQ7XG5cbmZ1bmN0aW9uIFVuaW9uRmluZChjb3VudCkge1xuICB0aGlzLnJvb3RzID0gbmV3IEFycmF5KGNvdW50KTtcbiAgdGhpcy5yYW5rcyA9IG5ldyBBcnJheShjb3VudCk7XG4gIFxuICBmb3IodmFyIGk9MDsgaTxjb3VudDsgKytpKSB7XG4gICAgdGhpcy5yb290c1tpXSA9IGk7XG4gICAgdGhpcy5yYW5rc1tpXSA9IDA7XG4gIH1cbn1cblxudmFyIHByb3RvID0gVW5pb25GaW5kLnByb3RvdHlwZVxuXG5PYmplY3QuZGVmaW5lUHJvcGVydHkocHJvdG8sIFwibGVuZ3RoXCIsIHtcbiAgXCJnZXRcIjogZnVuY3Rpb24oKSB7XG4gICAgcmV0dXJuIHRoaXMucm9vdHMubGVuZ3RoXG4gIH1cbn0pXG5cbnByb3RvLm1ha2VTZXQgPSBmdW5jdGlvbigpIHtcbiAgdmFyIG4gPSB0aGlzLnJvb3RzLmxlbmd0aDtcbiAgdGhpcy5yb290cy5wdXNoKG4pO1xuICB0aGlzLnJhbmtzLnB1c2goMCk7XG4gIHJldHVybiBuO1xufVxuXG5wcm90by5maW5kID0gZnVuY3Rpb24oeCkge1xuICB2YXIgcm9vdHMgPSB0aGlzLnJvb3RzO1xuICB3aGlsZShyb290c1t4XSAhPT0geCkge1xuICAgIHZhciB5ID0gcm9vdHNbeF07XG4gICAgcm9vdHNbeF0gPSByb290c1t5XTtcbiAgICB4ID0geTtcbiAgfVxuICByZXR1cm4geDtcbn1cblxucHJvdG8ubGluayA9IGZ1bmN0aW9uKHgsIHkpIHtcbiAgdmFyIHhyID0gdGhpcy5maW5kKHgpXG4gICAgLCB5ciA9IHRoaXMuZmluZCh5KTtcbiAgaWYoeHIgPT09IHlyKSB7XG4gICAgcmV0dXJuO1xuICB9XG4gIHZhciByYW5rcyA9IHRoaXMucmFua3NcbiAgICAsIHJvb3RzID0gdGhpcy5yb290c1xuICAgICwgeGQgICAgPSByYW5rc1t4cl1cbiAgICAsIHlkICAgID0gcmFua3NbeXJdO1xuICBpZih4ZCA8IHlkKSB7XG4gICAgcm9vdHNbeHJdID0geXI7XG4gIH0gZWxzZSBpZih5ZCA8IHhkKSB7XG4gICAgcm9vdHNbeXJdID0geHI7XG4gIH0gZWxzZSB7XG4gICAgcm9vdHNbeXJdID0geHI7XG4gICAgKytyYW5rc1t4cl07XG4gIH1cbn0iLCJcInVzZSBzdHJpY3RcIjsgXCJ1c2UgcmVzdHJpY3RcIjtcblxudmFyIGJpdHMgICAgICA9IHJlcXVpcmUoXCJiaXQtdHdpZGRsZVwiKVxuICAsIFVuaW9uRmluZCA9IHJlcXVpcmUoXCJ1bmlvbi1maW5kXCIpXG5cbi8vUmV0dXJucyB0aGUgZGltZW5zaW9uIG9mIGEgY2VsbCBjb21wbGV4XG5mdW5jdGlvbiBkaW1lbnNpb24oY2VsbHMpIHtcbiAgdmFyIGQgPSAwXG4gICAgLCBtYXggPSBNYXRoLm1heFxuICBmb3IodmFyIGk9MCwgaWw9Y2VsbHMubGVuZ3RoOyBpPGlsOyArK2kpIHtcbiAgICBkID0gbWF4KGQsIGNlbGxzW2ldLmxlbmd0aClcbiAgfVxuICByZXR1cm4gZC0xXG59XG5leHBvcnRzLmRpbWVuc2lvbiA9IGRpbWVuc2lvblxuXG4vL0NvdW50cyB0aGUgbnVtYmVyIG9mIHZlcnRpY2VzIGluIGZhY2VzXG5mdW5jdGlvbiBjb3VudFZlcnRpY2VzKGNlbGxzKSB7XG4gIHZhciB2YyA9IC0xXG4gICAgLCBtYXggPSBNYXRoLm1heFxuICBmb3IodmFyIGk9MCwgaWw9Y2VsbHMubGVuZ3RoOyBpPGlsOyArK2kpIHtcbiAgICB2YXIgYyA9IGNlbGxzW2ldXG4gICAgZm9yKHZhciBqPTAsIGpsPWMubGVuZ3RoOyBqPGpsOyArK2opIHtcbiAgICAgIHZjID0gbWF4KHZjLCBjW2pdKVxuICAgIH1cbiAgfVxuICByZXR1cm4gdmMrMVxufVxuZXhwb3J0cy5jb3VudFZlcnRpY2VzID0gY291bnRWZXJ0aWNlc1xuXG4vL1JldHVybnMgYSBkZWVwIGNvcHkgb2YgY2VsbHNcbmZ1bmN0aW9uIGNsb25lQ2VsbHMoY2VsbHMpIHtcbiAgdmFyIG5jZWxscyA9IG5ldyBBcnJheShjZWxscy5sZW5ndGgpXG4gIGZvcih2YXIgaT0wLCBpbD1jZWxscy5sZW5ndGg7IGk8aWw7ICsraSkge1xuICAgIG5jZWxsc1tpXSA9IGNlbGxzW2ldLnNsaWNlKDApXG4gIH1cbiAgcmV0dXJuIG5jZWxsc1xufVxuZXhwb3J0cy5jbG9uZUNlbGxzID0gY2xvbmVDZWxsc1xuXG4vL1JhbmtzIGEgcGFpciBvZiBjZWxscyB1cCB0byBwZXJtdXRhdGlvblxuZnVuY3Rpb24gY29tcGFyZUNlbGxzKGEsIGIpIHtcbiAgdmFyIG4gPSBhLmxlbmd0aFxuICAgICwgdCA9IGEubGVuZ3RoIC0gYi5sZW5ndGhcbiAgICAsIG1pbiA9IE1hdGgubWluXG4gIGlmKHQpIHtcbiAgICByZXR1cm4gdFxuICB9XG4gIHN3aXRjaChuKSB7XG4gICAgY2FzZSAwOlxuICAgICAgcmV0dXJuIDA7XG4gICAgY2FzZSAxOlxuICAgICAgcmV0dXJuIGFbMF0gLSBiWzBdO1xuICAgIGNhc2UgMjpcbiAgICAgIHZhciBkID0gYVswXSthWzFdLWJbMF0tYlsxXVxuICAgICAgaWYoZCkge1xuICAgICAgICByZXR1cm4gZFxuICAgICAgfVxuICAgICAgcmV0dXJuIG1pbihhWzBdLGFbMV0pIC0gbWluKGJbMF0sYlsxXSlcbiAgICBjYXNlIDM6XG4gICAgICB2YXIgbDEgPSBhWzBdK2FbMV1cbiAgICAgICAgLCBtMSA9IGJbMF0rYlsxXVxuICAgICAgZCA9IGwxK2FbMl0gLSAobTErYlsyXSlcbiAgICAgIGlmKGQpIHtcbiAgICAgICAgcmV0dXJuIGRcbiAgICAgIH1cbiAgICAgIHZhciBsMCA9IG1pbihhWzBdLCBhWzFdKVxuICAgICAgICAsIG0wID0gbWluKGJbMF0sIGJbMV0pXG4gICAgICAgICwgZCAgPSBtaW4obDAsIGFbMl0pIC0gbWluKG0wLCBiWzJdKVxuICAgICAgaWYoZCkge1xuICAgICAgICByZXR1cm4gZFxuICAgICAgfVxuICAgICAgcmV0dXJuIG1pbihsMCthWzJdLCBsMSkgLSBtaW4obTArYlsyXSwgbTEpXG4gICAgXG4gICAgLy9UT0RPOiBNYXliZSBvcHRpbWl6ZSBuPTQgYXMgd2VsbD9cbiAgICBcbiAgICBkZWZhdWx0OlxuICAgICAgdmFyIGFzID0gYS5zbGljZSgwKVxuICAgICAgYXMuc29ydCgpXG4gICAgICB2YXIgYnMgPSBiLnNsaWNlKDApXG4gICAgICBicy5zb3J0KClcbiAgICAgIGZvcih2YXIgaT0wOyBpPG47ICsraSkge1xuICAgICAgICB0ID0gYXNbaV0gLSBic1tpXVxuICAgICAgICBpZih0KSB7XG4gICAgICAgICAgcmV0dXJuIHRcbiAgICAgICAgfVxuICAgICAgfVxuICAgICAgcmV0dXJuIDBcbiAgfVxufVxuZXhwb3J0cy5jb21wYXJlQ2VsbHMgPSBjb21wYXJlQ2VsbHNcblxuZnVuY3Rpb24gY29tcGFyZVppcHBlZChhLCBiKSB7XG4gIHJldHVybiBjb21wYXJlQ2VsbHMoYVswXSwgYlswXSlcbn1cblxuLy9QdXRzIGEgY2VsbCBjb21wbGV4IGludG8gbm9ybWFsIG9yZGVyIGZvciB0aGUgcHVycG9zZXMgb2YgZmluZENlbGwgcXVlcmllc1xuZnVuY3Rpb24gbm9ybWFsaXplKGNlbGxzLCBhdHRyKSB7XG4gIGlmKGF0dHIpIHtcbiAgICB2YXIgbGVuID0gY2VsbHMubGVuZ3RoXG4gICAgdmFyIHppcHBlZCA9IG5ldyBBcnJheShsZW4pXG4gICAgZm9yKHZhciBpPTA7IGk8bGVuOyArK2kpIHtcbiAgICAgIHppcHBlZFtpXSA9IFtjZWxsc1tpXSwgYXR0cltpXV1cbiAgICB9XG4gICAgemlwcGVkLnNvcnQoY29tcGFyZVppcHBlZClcbiAgICBmb3IodmFyIGk9MDsgaTxsZW47ICsraSkge1xuICAgICAgY2VsbHNbaV0gPSB6aXBwZWRbaV1bMF1cbiAgICAgIGF0dHJbaV0gPSB6aXBwZWRbaV1bMV1cbiAgICB9XG4gICAgcmV0dXJuIGNlbGxzXG4gIH0gZWxzZSB7XG4gICAgY2VsbHMuc29ydChjb21wYXJlQ2VsbHMpXG4gICAgcmV0dXJuIGNlbGxzXG4gIH1cbn1cbmV4cG9ydHMubm9ybWFsaXplID0gbm9ybWFsaXplXG5cbi8vUmVtb3ZlcyBhbGwgZHVwbGljYXRlIGNlbGxzIGluIHRoZSBjb21wbGV4XG5mdW5jdGlvbiB1bmlxdWUoY2VsbHMpIHtcbiAgaWYoY2VsbHMubGVuZ3RoID09PSAwKSB7XG4gICAgcmV0dXJuIFtdXG4gIH1cbiAgdmFyIHB0ciA9IDFcbiAgICAsIGxlbiA9IGNlbGxzLmxlbmd0aFxuICBmb3IodmFyIGk9MTsgaTxsZW47ICsraSkge1xuICAgIHZhciBhID0gY2VsbHNbaV1cbiAgICBpZihjb21wYXJlQ2VsbHMoYSwgY2VsbHNbaS0xXSkpIHtcbiAgICAgIGlmKGkgPT09IHB0cikge1xuICAgICAgICBwdHIrK1xuICAgICAgICBjb250aW51ZVxuICAgICAgfVxuICAgICAgY2VsbHNbcHRyKytdID0gYVxuICAgIH1cbiAgfVxuICBjZWxscy5sZW5ndGggPSBwdHJcbiAgcmV0dXJuIGNlbGxzXG59XG5leHBvcnRzLnVuaXF1ZSA9IHVuaXF1ZTtcblxuLy9GaW5kcyBhIGNlbGwgaW4gYSBub3JtYWxpemVkIGNlbGwgY29tcGxleFxuZnVuY3Rpb24gZmluZENlbGwoY2VsbHMsIGMpIHtcbiAgdmFyIGxvID0gMFxuICAgICwgaGkgPSBjZWxscy5sZW5ndGgtMVxuICAgICwgciAgPSAtMVxuICB3aGlsZSAobG8gPD0gaGkpIHtcbiAgICB2YXIgbWlkID0gKGxvICsgaGkpID4+IDFcbiAgICAgICwgcyAgID0gY29tcGFyZUNlbGxzKGNlbGxzW21pZF0sIGMpXG4gICAgaWYocyA8PSAwKSB7XG4gICAgICBpZihzID09PSAwKSB7XG4gICAgICAgIHIgPSBtaWRcbiAgICAgIH1cbiAgICAgIGxvID0gbWlkICsgMVxuICAgIH0gZWxzZSBpZihzID4gMCkge1xuICAgICAgaGkgPSBtaWQgLSAxXG4gICAgfVxuICB9XG4gIHJldHVybiByXG59XG5leHBvcnRzLmZpbmRDZWxsID0gZmluZENlbGw7XG5cbi8vQnVpbGRzIGFuIGluZGV4IGZvciBhbiBuLWNlbGwuICBUaGlzIGlzIG1vcmUgZ2VuZXJhbCB0aGFuIGR1YWwsIGJ1dCBsZXNzIGVmZmljaWVudFxuZnVuY3Rpb24gaW5jaWRlbmNlKGZyb21fY2VsbHMsIHRvX2NlbGxzKSB7XG4gIHZhciBpbmRleCA9IG5ldyBBcnJheShmcm9tX2NlbGxzLmxlbmd0aClcbiAgZm9yKHZhciBpPTAsIGlsPWluZGV4Lmxlbmd0aDsgaTxpbDsgKytpKSB7XG4gICAgaW5kZXhbaV0gPSBbXVxuICB9XG4gIHZhciBiID0gW11cbiAgZm9yKHZhciBpPTAsIG49dG9fY2VsbHMubGVuZ3RoOyBpPG47ICsraSkge1xuICAgIHZhciBjID0gdG9fY2VsbHNbaV1cbiAgICB2YXIgY2wgPSBjLmxlbmd0aFxuICAgIGZvcih2YXIgaz0xLCBrbj0oMTw8Y2wpOyBrPGtuOyArK2spIHtcbiAgICAgIGIubGVuZ3RoID0gYml0cy5wb3BDb3VudChrKVxuICAgICAgdmFyIGwgPSAwXG4gICAgICBmb3IodmFyIGo9MDsgajxjbDsgKytqKSB7XG4gICAgICAgIGlmKGsgJiAoMTw8aikpIHtcbiAgICAgICAgICBiW2wrK10gPSBjW2pdXG4gICAgICAgIH1cbiAgICAgIH1cbiAgICAgIHZhciBpZHg9ZmluZENlbGwoZnJvbV9jZWxscywgYilcbiAgICAgIGlmKGlkeCA8IDApIHtcbiAgICAgICAgY29udGludWVcbiAgICAgIH1cbiAgICAgIHdoaWxlKHRydWUpIHtcbiAgICAgICAgaW5kZXhbaWR4KytdLnB1c2goaSlcbiAgICAgICAgaWYoaWR4ID49IGZyb21fY2VsbHMubGVuZ3RoIHx8IGNvbXBhcmVDZWxscyhmcm9tX2NlbGxzW2lkeF0sIGIpICE9PSAwKSB7XG4gICAgICAgICAgYnJlYWtcbiAgICAgICAgfVxuICAgICAgfVxuICAgIH1cbiAgfVxuICByZXR1cm4gaW5kZXhcbn1cbmV4cG9ydHMuaW5jaWRlbmNlID0gaW5jaWRlbmNlXG5cbi8vQ29tcHV0ZXMgdGhlIGR1YWwgb2YgdGhlIG1lc2guICBUaGlzIGlzIGJhc2ljYWxseSBhbiBvcHRpbWl6ZWQgdmVyc2lvbiBvZiBidWlsZEluZGV4IGZvciB0aGUgc2l0dWF0aW9uIHdoZXJlIGZyb21fY2VsbHMgaXMganVzdCB0aGUgbGlzdCBvZiB2ZXJ0aWNlc1xuZnVuY3Rpb24gZHVhbChjZWxscywgdmVydGV4X2NvdW50KSB7XG4gIGlmKCF2ZXJ0ZXhfY291bnQpIHtcbiAgICByZXR1cm4gaW5jaWRlbmNlKHVuaXF1ZShza2VsZXRvbihjZWxscywgMCkpLCBjZWxscywgMClcbiAgfVxuICB2YXIgcmVzID0gbmV3IEFycmF5KHZlcnRleF9jb3VudClcbiAgZm9yKHZhciBpPTA7IGk8dmVydGV4X2NvdW50OyArK2kpIHtcbiAgICByZXNbaV0gPSBbXVxuICB9XG4gIGZvcih2YXIgaT0wLCBsZW49Y2VsbHMubGVuZ3RoOyBpPGxlbjsgKytpKSB7XG4gICAgdmFyIGMgPSBjZWxsc1tpXVxuICAgIGZvcih2YXIgaj0wLCBjbD1jLmxlbmd0aDsgajxjbDsgKytqKSB7XG4gICAgICByZXNbY1tqXV0ucHVzaChpKVxuICAgIH1cbiAgfVxuICByZXR1cm4gcmVzXG59XG5leHBvcnRzLmR1YWwgPSBkdWFsXG5cbi8vRW51bWVyYXRlcyBhbGwgY2VsbHMgaW4gdGhlIGNvbXBsZXhcbmZ1bmN0aW9uIGV4cGxvZGUoY2VsbHMpIHtcbiAgdmFyIHJlc3VsdCA9IFtdXG4gIGZvcih2YXIgaT0wLCBpbD1jZWxscy5sZW5ndGg7IGk8aWw7ICsraSkge1xuICAgIHZhciBjID0gY2VsbHNbaV1cbiAgICAgICwgY2wgPSBjLmxlbmd0aHwwXG4gICAgZm9yKHZhciBqPTEsIGpsPSgxPDxjbCk7IGo8amw7ICsraikge1xuICAgICAgdmFyIGIgPSBbXVxuICAgICAgZm9yKHZhciBrPTA7IGs8Y2w7ICsraykge1xuICAgICAgICBpZigoaiA+Pj4gaykgJiAxKSB7XG4gICAgICAgICAgYi5wdXNoKGNba10pXG4gICAgICAgIH1cbiAgICAgIH1cbiAgICAgIHJlc3VsdC5wdXNoKGIpXG4gICAgfVxuICB9XG4gIHJldHVybiBub3JtYWxpemUocmVzdWx0KVxufVxuZXhwb3J0cy5leHBsb2RlID0gZXhwbG9kZVxuXG4vL0VudW1lcmF0ZXMgYWxsIG9mIHRoZSBuLWNlbGxzIG9mIGEgY2VsbCBjb21wbGV4XG5mdW5jdGlvbiBza2VsZXRvbihjZWxscywgbikge1xuICBpZihuIDwgMCkge1xuICAgIHJldHVybiBbXVxuICB9XG4gIHZhciByZXN1bHQgPSBbXVxuICAgICwgazAgICAgID0gKDE8PChuKzEpKS0xXG4gIGZvcih2YXIgaT0wOyBpPGNlbGxzLmxlbmd0aDsgKytpKSB7XG4gICAgdmFyIGMgPSBjZWxsc1tpXVxuICAgIGZvcih2YXIgaz1rMDsgazwoMTw8Yy5sZW5ndGgpOyBrPWJpdHMubmV4dENvbWJpbmF0aW9uKGspKSB7XG4gICAgICB2YXIgYiA9IG5ldyBBcnJheShuKzEpXG4gICAgICAgICwgbCA9IDBcbiAgICAgIGZvcih2YXIgaj0wOyBqPGMubGVuZ3RoOyArK2opIHtcbiAgICAgICAgaWYoayAmICgxPDxqKSkge1xuICAgICAgICAgIGJbbCsrXSA9IGNbal1cbiAgICAgICAgfVxuICAgICAgfVxuICAgICAgcmVzdWx0LnB1c2goYilcbiAgICB9XG4gIH1cbiAgcmV0dXJuIG5vcm1hbGl6ZShyZXN1bHQpXG59XG5leHBvcnRzLnNrZWxldG9uID0gc2tlbGV0b247XG5cbi8vQ29tcHV0ZXMgdGhlIGJvdW5kYXJ5IG9mIGFsbCBjZWxscywgZG9lcyBub3QgcmVtb3ZlIGR1cGxpY2F0ZXNcbmZ1bmN0aW9uIGJvdW5kYXJ5KGNlbGxzKSB7XG4gIHZhciByZXMgPSBbXVxuICBmb3IodmFyIGk9MCxpbD1jZWxscy5sZW5ndGg7IGk8aWw7ICsraSkge1xuICAgIHZhciBjID0gY2VsbHNbaV1cbiAgICBmb3IodmFyIGo9MCxjbD1jLmxlbmd0aDsgajxjbDsgKytqKSB7XG4gICAgICB2YXIgYiA9IG5ldyBBcnJheShjLmxlbmd0aC0xKVxuICAgICAgZm9yKHZhciBrPTAsIGw9MDsgazxjbDsgKytrKSB7XG4gICAgICAgIGlmKGsgIT09IGopIHtcbiAgICAgICAgICBiW2wrK10gPSBjW2tdXG4gICAgICAgIH1cbiAgICAgIH1cbiAgICAgIHJlcy5wdXNoKGIpXG4gICAgfVxuICB9XG4gIHJldHVybiBub3JtYWxpemUocmVzKVxufVxuZXhwb3J0cy5ib3VuZGFyeSA9IGJvdW5kYXJ5O1xuXG4vL0NvbXB1dGVzIGNvbm5lY3RlZCBjb21wb25lbnRzIGZvciBhIGRlbnNlIGNlbGwgY29tcGxleFxuZnVuY3Rpb24gY29ubmVjdGVkQ29tcG9uZW50c19kZW5zZShjZWxscywgdmVydGV4X2NvdW50KSB7XG4gIHZhciBsYWJlbHMgPSBuZXcgVW5pb25GaW5kKHZlcnRleF9jb3VudClcbiAgZm9yKHZhciBpPTA7IGk8Y2VsbHMubGVuZ3RoOyArK2kpIHtcbiAgICB2YXIgYyA9IGNlbGxzW2ldXG4gICAgZm9yKHZhciBqPTA7IGo8Yy5sZW5ndGg7ICsraikge1xuICAgICAgZm9yKHZhciBrPWorMTsgazxjLmxlbmd0aDsgKytrKSB7XG4gICAgICAgIGxhYmVscy5saW5rKGNbal0sIGNba10pXG4gICAgICB9XG4gICAgfVxuICB9XG4gIHZhciBjb21wb25lbnRzID0gW11cbiAgICAsIGNvbXBvbmVudF9sYWJlbHMgPSBsYWJlbHMucmFua3NcbiAgZm9yKHZhciBpPTA7IGk8Y29tcG9uZW50X2xhYmVscy5sZW5ndGg7ICsraSkge1xuICAgIGNvbXBvbmVudF9sYWJlbHNbaV0gPSAtMVxuICB9XG4gIGZvcih2YXIgaT0wOyBpPGNlbGxzLmxlbmd0aDsgKytpKSB7XG4gICAgdmFyIGwgPSBsYWJlbHMuZmluZChjZWxsc1tpXVswXSlcbiAgICBpZihjb21wb25lbnRfbGFiZWxzW2xdIDwgMCkge1xuICAgICAgY29tcG9uZW50X2xhYmVsc1tsXSA9IGNvbXBvbmVudHMubGVuZ3RoXG4gICAgICBjb21wb25lbnRzLnB1c2goW2NlbGxzW2ldLnNsaWNlKDApXSlcbiAgICB9IGVsc2Uge1xuICAgICAgY29tcG9uZW50c1tjb21wb25lbnRfbGFiZWxzW2xdXS5wdXNoKGNlbGxzW2ldLnNsaWNlKDApKVxuICAgIH1cbiAgfVxuICByZXR1cm4gY29tcG9uZW50c1xufVxuXG4vL0NvbXB1dGVzIGNvbm5lY3RlZCBjb21wb25lbnRzIGZvciBhIHNwYXJzZSBncmFwaFxuZnVuY3Rpb24gY29ubmVjdGVkQ29tcG9uZW50c19zcGFyc2UoY2VsbHMpIHtcbiAgdmFyIHZlcnRpY2VzICA9IHVuaXF1ZShub3JtYWxpemUoc2tlbGV0b24oY2VsbHMsIDApKSlcbiAgICAsIGxhYmVscyAgICA9IG5ldyBVbmlvbkZpbmQodmVydGljZXMubGVuZ3RoKVxuICBmb3IodmFyIGk9MDsgaTxjZWxscy5sZW5ndGg7ICsraSkge1xuICAgIHZhciBjID0gY2VsbHNbaV1cbiAgICBmb3IodmFyIGo9MDsgajxjLmxlbmd0aDsgKytqKSB7XG4gICAgICB2YXIgdmogPSBmaW5kQ2VsbCh2ZXJ0aWNlcywgW2Nbal1dKVxuICAgICAgZm9yKHZhciBrPWorMTsgazxjLmxlbmd0aDsgKytrKSB7XG4gICAgICAgIGxhYmVscy5saW5rKHZqLCBmaW5kQ2VsbCh2ZXJ0aWNlcywgW2Nba11dKSlcbiAgICAgIH1cbiAgICB9XG4gIH1cbiAgdmFyIGNvbXBvbmVudHMgICAgICAgID0gW11cbiAgICAsIGNvbXBvbmVudF9sYWJlbHMgID0gbGFiZWxzLnJhbmtzXG4gIGZvcih2YXIgaT0wOyBpPGNvbXBvbmVudF9sYWJlbHMubGVuZ3RoOyArK2kpIHtcbiAgICBjb21wb25lbnRfbGFiZWxzW2ldID0gLTFcbiAgfVxuICBmb3IodmFyIGk9MDsgaTxjZWxscy5sZW5ndGg7ICsraSkge1xuICAgIHZhciBsID0gbGFiZWxzLmZpbmQoZmluZENlbGwodmVydGljZXMsIFtjZWxsc1tpXVswXV0pKTtcbiAgICBpZihjb21wb25lbnRfbGFiZWxzW2xdIDwgMCkge1xuICAgICAgY29tcG9uZW50X2xhYmVsc1tsXSA9IGNvbXBvbmVudHMubGVuZ3RoXG4gICAgICBjb21wb25lbnRzLnB1c2goW2NlbGxzW2ldLnNsaWNlKDApXSlcbiAgICB9IGVsc2Uge1xuICAgICAgY29tcG9uZW50c1tjb21wb25lbnRfbGFiZWxzW2xdXS5wdXNoKGNlbGxzW2ldLnNsaWNlKDApKVxuICAgIH1cbiAgfVxuICByZXR1cm4gY29tcG9uZW50c1xufVxuXG4vL0NvbXB1dGVzIGNvbm5lY3RlZCBjb21wb25lbnRzIGZvciBhIGNlbGwgY29tcGxleFxuZnVuY3Rpb24gY29ubmVjdGVkQ29tcG9uZW50cyhjZWxscywgdmVydGV4X2NvdW50KSB7XG4gIGlmKHZlcnRleF9jb3VudCkge1xuICAgIHJldHVybiBjb25uZWN0ZWRDb21wb25lbnRzX2RlbnNlKGNlbGxzLCB2ZXJ0ZXhfY291bnQpXG4gIH1cbiAgcmV0dXJuIGNvbm5lY3RlZENvbXBvbmVudHNfc3BhcnNlKGNlbGxzKVxufVxuZXhwb3J0cy5jb25uZWN0ZWRDb21wb25lbnRzID0gY29ubmVjdGVkQ29tcG9uZW50c1xuIiwiXCJ1c2Ugc3RyaWN0XCJcblxuZnVuY3Rpb24gdW5pcXVlX3ByZWQobGlzdCwgY29tcGFyZSkge1xuICB2YXIgcHRyID0gMVxuICAgICwgbGVuID0gbGlzdC5sZW5ndGhcbiAgICAsIGE9bGlzdFswXSwgYj1saXN0WzBdXG4gIGZvcih2YXIgaT0xOyBpPGxlbjsgKytpKSB7XG4gICAgYiA9IGFcbiAgICBhID0gbGlzdFtpXVxuICAgIGlmKGNvbXBhcmUoYSwgYikpIHtcbiAgICAgIGlmKGkgPT09IHB0cikge1xuICAgICAgICBwdHIrK1xuICAgICAgICBjb250aW51ZVxuICAgICAgfVxuICAgICAgbGlzdFtwdHIrK10gPSBhXG4gICAgfVxuICB9XG4gIGxpc3QubGVuZ3RoID0gcHRyXG4gIHJldHVybiBsaXN0XG59XG5cbmZ1bmN0aW9uIHVuaXF1ZV9lcShsaXN0KSB7XG4gIHZhciBwdHIgPSAxXG4gICAgLCBsZW4gPSBsaXN0Lmxlbmd0aFxuICAgICwgYT1saXN0WzBdLCBiID0gbGlzdFswXVxuICBmb3IodmFyIGk9MTsgaTxsZW47ICsraSwgYj1hKSB7XG4gICAgYiA9IGFcbiAgICBhID0gbGlzdFtpXVxuICAgIGlmKGEgIT09IGIpIHtcbiAgICAgIGlmKGkgPT09IHB0cikge1xuICAgICAgICBwdHIrK1xuICAgICAgICBjb250aW51ZVxuICAgICAgfVxuICAgICAgbGlzdFtwdHIrK10gPSBhXG4gICAgfVxuICB9XG4gIGxpc3QubGVuZ3RoID0gcHRyXG4gIHJldHVybiBsaXN0XG59XG5cbmZ1bmN0aW9uIHVuaXF1ZShsaXN0LCBjb21wYXJlLCBzb3J0ZWQpIHtcbiAgaWYobGlzdC5sZW5ndGggPT09IDApIHtcbiAgICByZXR1cm4gbGlzdFxuICB9XG4gIGlmKGNvbXBhcmUpIHtcbiAgICBpZighc29ydGVkKSB7XG4gICAgICBsaXN0LnNvcnQoY29tcGFyZSlcbiAgICB9XG4gICAgcmV0dXJuIHVuaXF1ZV9wcmVkKGxpc3QsIGNvbXBhcmUpXG4gIH1cbiAgaWYoIXNvcnRlZCkge1xuICAgIGxpc3Quc29ydCgpXG4gIH1cbiAgcmV0dXJuIHVuaXF1ZV9lcShsaXN0KVxufVxuXG5tb2R1bGUuZXhwb3J0cyA9IHVuaXF1ZVxuIiwiXCJ1c2Ugc3RyaWN0XCJcblxudmFyIGNoID0gcmVxdWlyZShcImluY3JlbWVudGFsLWNvbnZleC1odWxsXCIpXG52YXIgdW5pcSA9IHJlcXVpcmUoXCJ1bmlxXCIpXG5cbm1vZHVsZS5leHBvcnRzID0gdHJpYW5ndWxhdGVcblxuZnVuY3Rpb24gTGlmdGVkUG9pbnQocCwgaSkge1xuICB0aGlzLnBvaW50ID0gcFxuICB0aGlzLmluZGV4ID0gaVxufVxuXG5mdW5jdGlvbiBjb21wYXJlTGlmdGVkKGEsIGIpIHtcbiAgdmFyIGFwID0gYS5wb2ludFxuICB2YXIgYnAgPSBiLnBvaW50XG4gIHZhciBkID0gYXAubGVuZ3RoXG4gIGZvcih2YXIgaT0wOyBpPGQ7ICsraSkge1xuICAgIHZhciBzID0gYnBbaV0gLSBhcFtpXVxuICAgIGlmKHMpIHtcbiAgICAgIHJldHVybiBzXG4gICAgfVxuICB9XG4gIHJldHVybiAwXG59XG5cbmZ1bmN0aW9uIHRyaWFuZ3VsYXRlMUQobiwgcG9pbnRzLCBpbmNsdWRlUG9pbnRBdEluZmluaXR5KSB7XG4gIGlmKG4gPT09IDEpIHtcbiAgICBpZihpbmNsdWRlUG9pbnRBdEluZmluaXR5KSB7XG4gICAgICByZXR1cm4gWyBbLTEsIDBdIF1cbiAgICB9IGVsc2Uge1xuICAgICAgcmV0dXJuIFtdXG4gICAgfVxuICB9XG4gIHZhciBsaWZ0ZWQgPSBwb2ludHMubWFwKGZ1bmN0aW9uKHAsIGkpIHtcbiAgICByZXR1cm4gWyBwWzBdLCBpIF1cbiAgfSlcbiAgbGlmdGVkLnNvcnQoZnVuY3Rpb24oYSxiKSB7XG4gICAgcmV0dXJuIGFbMF0gLSBiWzBdXG4gIH0pXG4gIHZhciBjZWxscyA9IG5ldyBBcnJheShuIC0gMSlcbiAgZm9yKHZhciBpPTE7IGk8bjsgKytpKSB7XG4gICAgdmFyIGEgPSBsaWZ0ZWRbaS0xXVxuICAgIHZhciBiID0gbGlmdGVkW2ldXG4gICAgY2VsbHNbaS0xXSA9IFsgYVsxXSwgYlsxXSBdXG4gIH1cbiAgaWYoaW5jbHVkZVBvaW50QXRJbmZpbml0eSkge1xuICAgIGNlbGxzLnB1c2goXG4gICAgICBbIC0xLCBjZWxsc1swXVsxXSwgXSxcbiAgICAgIFsgY2VsbHNbbi0xXVsxXSwgLTEgXSlcbiAgfVxuICByZXR1cm4gY2VsbHNcbn1cblxuZnVuY3Rpb24gdHJpYW5ndWxhdGUocG9pbnRzLCBpbmNsdWRlUG9pbnRBdEluZmluaXR5KSB7XG4gIHZhciBuID0gcG9pbnRzLmxlbmd0aFxuICBpZihuID09PSAwKSB7XG4gICAgcmV0dXJuIFtdXG4gIH1cbiAgXG4gIHZhciBkID0gcG9pbnRzWzBdLmxlbmd0aFxuICBpZihkIDwgMSkge1xuICAgIHJldHVybiBbXVxuICB9XG5cbiAgLy9TcGVjaWFsIGNhc2U6ICBGb3IgMUQgd2UgY2FuIGp1c3Qgc29ydCB0aGUgcG9pbnRzXG4gIGlmKGQgPT09IDEpIHtcbiAgICByZXR1cm4gdHJpYW5ndWxhdGUxRChuLCBwb2ludHMsIGluY2x1ZGVQb2ludEF0SW5maW5pdHkpXG4gIH1cbiAgXG4gIC8vTGlmdCBwb2ludHMsIHNvcnRcbiAgdmFyIGxpZnRlZCA9IG5ldyBBcnJheShuKVxuICB2YXIgdXBwZXIgPSAxLjBcbiAgZm9yKHZhciBpPTA7IGk8bjsgKytpKSB7XG4gICAgdmFyIHAgPSBwb2ludHNbaV1cbiAgICB2YXIgeCA9IG5ldyBBcnJheShkKzEpXG4gICAgdmFyIGwgPSAwLjBcbiAgICBmb3IodmFyIGo9MDsgajxkOyArK2opIHtcbiAgICAgIHZhciB2ID0gcFtqXVxuICAgICAgeFtqXSA9IHZcbiAgICAgIGwgKz0gdiAqIHZcbiAgICB9XG4gICAgeFtkXSA9IGxcbiAgICBsaWZ0ZWRbaV0gPSBuZXcgTGlmdGVkUG9pbnQoeCwgaSlcbiAgICB1cHBlciA9IE1hdGgubWF4KGwsIHVwcGVyKVxuICB9XG4gIHVuaXEobGlmdGVkLCBjb21wYXJlTGlmdGVkKVxuICBcbiAgLy9Eb3VibGUgcG9pbnRzXG4gIG4gPSBsaWZ0ZWQubGVuZ3RoXG5cbiAgLy9DcmVhdGUgbmV3IGxpc3Qgb2YgcG9pbnRzXG4gIHZhciBkcG9pbnRzID0gbmV3IEFycmF5KG4gKyBkICsgMSlcbiAgdmFyIGRpbmRleCA9IG5ldyBBcnJheShuICsgZCArIDEpXG5cbiAgLy9BZGQgc3RlaW5lciBwb2ludHMgYXQgdG9wXG4gIHZhciB1ID0gKGQrMSkgKiAoZCsxKSAqIHVwcGVyXG4gIHZhciB5ID0gbmV3IEFycmF5KGQrMSlcbiAgZm9yKHZhciBpPTA7IGk8PWQ7ICsraSkge1xuICAgIHlbaV0gPSAwLjBcbiAgfVxuICB5W2RdID0gdVxuXG4gIGRwb2ludHNbMF0gPSB5LnNsaWNlKClcbiAgZGluZGV4WzBdID0gLTFcblxuICBmb3IodmFyIGk9MDsgaTw9ZDsgKytpKSB7XG4gICAgdmFyIHggPSB5LnNsaWNlKClcbiAgICB4W2ldID0gMVxuICAgIGRwb2ludHNbaSsxXSA9IHhcbiAgICBkaW5kZXhbaSsxXSA9IC0xXG4gIH1cblxuICAvL0NvcHkgcmVzdCBvZiB0aGUgcG9pbnRzIG92ZXJcbiAgZm9yKHZhciBpPTA7IGk8bjsgKytpKSB7XG4gICAgdmFyIGggPSBsaWZ0ZWRbaV1cbiAgICBkcG9pbnRzW2kgKyBkICsgMV0gPSBoLnBvaW50XG4gICAgZGluZGV4W2kgKyBkICsgMV0gPSAgaC5pbmRleFxuICB9XG5cbiAgLy9Db25zdHJ1Y3QgY29udmV4IGh1bGxcbiAgdmFyIGh1bGwgPSBjaChkcG9pbnRzLCBmYWxzZSlcbiAgaWYoaW5jbHVkZVBvaW50QXRJbmZpbml0eSkge1xuICAgIGh1bGwgPSBodWxsLmZpbHRlcihmdW5jdGlvbihjZWxsKSB7XG4gICAgICB2YXIgY291bnQgPSAwXG4gICAgICBmb3IodmFyIGo9MDsgajw9ZDsgKytqKSB7XG4gICAgICAgIHZhciB2ID0gZGluZGV4W2NlbGxbal1dXG4gICAgICAgIGlmKHYgPCAwKSB7XG4gICAgICAgICAgaWYoKytjb3VudCA+PSAyKSB7XG4gICAgICAgICAgICByZXR1cm4gZmFsc2VcbiAgICAgICAgICB9XG4gICAgICAgIH1cbiAgICAgICAgY2VsbFtqXSA9IHZcbiAgICAgIH1cbiAgICAgIHJldHVybiB0cnVlXG4gICAgfSlcbiAgfSBlbHNlIHtcbiAgICBodWxsID0gaHVsbC5maWx0ZXIoZnVuY3Rpb24oY2VsbCkge1xuICAgICAgZm9yKHZhciBpPTA7IGk8PWQ7ICsraSkge1xuICAgICAgICB2YXIgdiA9IGRpbmRleFtjZWxsW2ldXVxuICAgICAgICBpZih2IDwgMCkge1xuICAgICAgICAgIHJldHVybiBmYWxzZVxuICAgICAgICB9XG4gICAgICAgIGNlbGxbaV0gPSB2XG4gICAgICB9XG4gICAgICByZXR1cm4gdHJ1ZVxuICAgIH0pXG4gIH1cblxuICBpZihkICYgMSkge1xuICAgIGZvcih2YXIgaT0wOyBpPGh1bGwubGVuZ3RoOyArK2kpIHtcbiAgICAgIHZhciBoID0gaHVsbFtpXVxuICAgICAgdmFyIHggPSBoWzBdXG4gICAgICBoWzBdID0gaFsxXVxuICAgICAgaFsxXSA9IHhcbiAgICB9XG4gIH1cblxuICByZXR1cm4gaHVsbFxufSIsIlxubW9kdWxlLmV4cG9ydHMgPSBwYXJzZVxuXG4vKipcbiAqIGV4cGVjdGVkIGFyZ3VtZW50IGxlbmd0aHNcbiAqIEB0eXBlIHtPYmplY3R9XG4gKi9cblxudmFyIGxlbmd0aCA9IHthOiA3LCBjOiA2LCBoOiAxLCBsOiAyLCBtOiAyLCBxOiA0LCBzOiA0LCB0OiAyLCB2OiAxLCB6OiAwfVxuXG4vKipcbiAqIHNlZ21lbnQgcGF0dGVyblxuICogQHR5cGUge1JlZ0V4cH1cbiAqL1xuXG52YXIgc2VnbWVudCA9IC8oW2FzdHZ6cW1obGNdKShbXmFzdHZ6cW1obGNdKikvaWdcblxuLyoqXG4gKiBwYXJzZSBhbiBzdmcgcGF0aCBkYXRhIHN0cmluZy4gR2VuZXJhdGVzIGFuIEFycmF5XG4gKiBvZiBjb21tYW5kcyB3aGVyZSBlYWNoIGNvbW1hbmQgaXMgYW4gQXJyYXkgb2YgdGhlXG4gKiBmb3JtIGBbY29tbWFuZCwgYXJnMSwgYXJnMiwgLi4uXWBcbiAqXG4gKiBAcGFyYW0ge1N0cmluZ30gcGF0aFxuICogQHJldHVybiB7QXJyYXl9XG4gKi9cblxuZnVuY3Rpb24gcGFyc2UocGF0aCkge1xuXHR2YXIgZGF0YSA9IFtdXG5cdHBhdGgucmVwbGFjZShzZWdtZW50LCBmdW5jdGlvbihfLCBjb21tYW5kLCBhcmdzKXtcblx0XHR2YXIgdHlwZSA9IGNvbW1hbmQudG9Mb3dlckNhc2UoKVxuXHRcdGFyZ3MgPSBwYXJzZVZhbHVlcyhhcmdzKVxuXG5cdFx0Ly8gb3ZlcmxvYWRlZCBtb3ZlVG9cblx0XHRpZiAodHlwZSA9PSAnbScgJiYgYXJncy5sZW5ndGggPiAyKSB7XG5cdFx0XHRkYXRhLnB1c2goW2NvbW1hbmRdLmNvbmNhdChhcmdzLnNwbGljZSgwLCAyKSkpXG5cdFx0XHR0eXBlID0gJ2wnXG5cdFx0XHRjb21tYW5kID0gY29tbWFuZCA9PSAnbScgPyAnbCcgOiAnTCdcblx0XHR9XG5cblx0XHR3aGlsZSAodHJ1ZSkge1xuXHRcdFx0aWYgKGFyZ3MubGVuZ3RoID09IGxlbmd0aFt0eXBlXSkge1xuXHRcdFx0XHRhcmdzLnVuc2hpZnQoY29tbWFuZClcblx0XHRcdFx0cmV0dXJuIGRhdGEucHVzaChhcmdzKVxuXHRcdFx0fVxuXHRcdFx0aWYgKGFyZ3MubGVuZ3RoIDwgbGVuZ3RoW3R5cGVdKSB0aHJvdyBuZXcgRXJyb3IoJ21hbGZvcm1lZCBwYXRoIGRhdGEnKVxuXHRcdFx0ZGF0YS5wdXNoKFtjb21tYW5kXS5jb25jYXQoYXJncy5zcGxpY2UoMCwgbGVuZ3RoW3R5cGVdKSkpXG5cdFx0fVxuXHR9KVxuXHRyZXR1cm4gZGF0YVxufVxuXG5mdW5jdGlvbiBwYXJzZVZhbHVlcyhhcmdzKXtcblx0YXJncyA9IGFyZ3MubWF0Y2goLy0/Wy4wLTldKyg/OmVbLStdP1xcZCspPy9pZylcblx0cmV0dXJuIGFyZ3MgPyBhcmdzLm1hcChOdW1iZXIpIDogW11cbn1cbiIsIid1c2Ugc3RyaWN0JztcblxudmFyIHNpZ24gPSByZXF1aXJlKCcuL3V0aWxpdGllcy5qcycpLnNpZ247XG52YXIgY2FsY3VsYXRlRGlzdGFuY2UgPSByZXF1aXJlKCcuL3V0aWxpdGllcy5qcycpLmRpc3RhbmNlO1xuXG52YXIgcG9pbnRzID0gcmVxdWlyZSgnLi9pbml0aWFsaXplUG9pbnRzLmpzJykucG9pbnRzO1xudmFyIGNpdHlTZXQgPSByZXF1aXJlKCcuL2luaXRpYWxpemVQb2ludHMuanMnKS5jaXR5U2V0O1xudmFyIHRleHRQb2ludHNJZCA9IHJlcXVpcmUoJy4vaW5pdGlhbGl6ZVBvaW50cy5qcycpLnRleHRQb2ludHNJZDtcbnZhciBwb3NzaWJsZVN0YXJ0UG9pbnRzSWQgPSByZXF1aXJlKCcuL2luaXRpYWxpemVQb2ludHMuanMnKS5wb3NzaWJsZVN0YXJ0UG9pbnRzSWQ7XG5cbnZhciBsaXZlTW91c2VQb3NpdGlvbiA9IHJlcXVpcmUoJy4vbW91c2UuanMnKTtcblxudmFyIFZlY3RvciA9IHJlcXVpcmUoJy4vdmVjdG9yLmpzJyk7XG5cbnZhciByYW5kb20gPSBNYXRoLnJhbmRvbTtcbnZhciBmbG9vciA9IE1hdGguZmxvb3I7XG52YXIgUkVQVUxTSU9OID0gMC4wNTtcblxuXG5tb2R1bGUuZXhwb3J0cyA9IGZ1bmN0aW9uKGNvbnRhaW5lcil7XG5cbiAgICB2YXIgbW91c2UgPSBsaXZlTW91c2VQb3NpdGlvbihjb250YWluZXIpO1xuXG5cbiAgICBmdW5jdGlvbiBBbnQocG9pbnQpIHtcbiAgICAgICAgdGhpcy54ID0gcG9pbnQueDsgICAgICAgICAgICAgICAgXG4gICAgICAgIHRoaXMueSA9IHBvaW50Lnk7XG4gICAgICAgIHRoaXMudmVsb2NpdHkgPSAwLjAwMztcbiAgICAgICAgdGhpcy5lZGdlID0gdW5kZWZpbmVkO1xuICAgICAgICB0aGlzLnN0YXRlID0gXCJmb3JhZ2VcIjtcbiAgICAgICAgdGhpcy5lZGdlcyA9IFtdO1xuICAgICAgICB0aGlzLmxhc3RDaXR5ID0gdW5kZWZpbmVkO1xuICAgICAgICB0aGlzLm9yaWdpbiA9IHBvaW50O1xuICAgICAgICB0aGlzLmRlc3RpbmF0aW9uID0gdW5kZWZpbmVkO1xuICAgICAgICB0aGlzLm9yaWVudGF0aW9uID0gdW5kZWZpbmVkO1xuICAgICAgICB0aGlzLmRpcmVjdGlvbiA9IG5ldyBWZWN0b3IoMCwwKTtcbiAgICAgICAgdGhpcy5wcm9nID0gMDtcbiAgICB9XG4gICAgLy8gZm9yYWdlOiB0aGUgYW50IHdhbmRlcnMgYXJvdW5kIHdpdGhvdXQgYW55IHBoZXJvbW9uIGRlcG9zaXRpb25cbiAgICAvLyBvbmNlIGl0IGZpbmRzIGEgY2l0eSwgaXQgc3RhcnRzIHJlbWVtYmVyaW5nIHRoZSBub2RlcyBpdCBnb2VzIHRocm91Z2hcbiAgICAvLyB3aGVuIGl0IGZpbmRzIGFub3RoZXIgY2l0eSwgaXQgY29tcHV0ZXMgdGhlIHBhdGggbGVuZ3RoIGFuZCBhZGRzIHBoZXJvbW9ucyBvbmUgZWFjaCBlZGdlc1xuICAgIC8vIHByb3BvcnRpb25uYWx5IHRvIHRoZSBzaG9ydGVzdG5lc3Mgb2YgdGhlIHBhdGhcbiAgICAvLyBpdCByZXNldHMgdGhlIGxpc3Qgb2Ygbm9kZXMgYW5kIGNvbnRpbnVlc1xuICAgIC8vIHdoaWxlIGZvcmFnaW5nIHRoZSBhbnQgY2hvc2VzIHRoZSBwYXRoIHdpdGggYSBwaGVyb21vbiBwcmVmZXJlbmNlXG5cblxuICAgIC8vIHN0YXRpYyBtZXRob2RzXG4gICAgQW50LmdlbmVyYXRlUmFuZFN0YXJ0UG9pbnQgPSBmdW5jdGlvbigpIHtcbiAgICAgICAgdmFyIHJhbmRJZCA9IE1hdGguZmxvb3IocG9zc2libGVTdGFydFBvaW50c0lkLmxlbmd0aCAqIHJhbmRvbSgpKTtcbiAgICAgICAgdmFyIHJhbmRTdGFydFBvaW50ID0gcG9pbnRzW3Bvc3NpYmxlU3RhcnRQb2ludHNJZFtyYW5kSWRdXTtcbiAgICAgICAgcmV0dXJuIHJhbmRTdGFydFBvaW50O1xuICAgIH1cblxuXG4gICAgLy8gbWV0aG9kc1xuICAgIEFudC5wcm90b3R5cGUgPSB7XG5cbiAgICAgICAgdHJhbnNpdDogZnVuY3Rpb24oKXtcbiAgICAgICAgICAgIHN3aXRjaCAodGhpcy5zdGF0ZSkge1xuICAgICAgICAgICAgY2FzZSBcImZvcmFnZVwiOlxuICAgICAgICAgICAgICAgIHZhciByZXMgPSB0aGlzLm1vdmUoKTtcbiAgICAgICAgICAgICAgICBpZiAocmVzLmNpdHlSZWFjaGVkKSB7XG4gICAgICAgICAgICAgICAgICAgIHRoaXMuc3RhdGUgPSBcInBoZXJvbW9uaW5nXCI7XG4gICAgICAgICAgICAgICAgICAgIHRoaXMubGFzdENpdHkgPSB0aGlzLm9yaWdpbi5pZDtcbiAgICAgICAgICAgICAgICB9O1xuICAgICAgICAgICAgICAgIGJyZWFrO1xuICAgICAgICAgICAgY2FzZSBcInBoZXJvbW9uaW5nXCI6XG4gICAgICAgICAgICAgICAgdmFyIHJlcyA9IHRoaXMubW92ZSgpO1xuICAgICAgICAgICAgICAgIGlmIChyZXMuZWRnZUNoYW5nZWQpIHtcbiAgICAgICAgICAgICAgICAgICAgdGhpcy5lZGdlcy5wdXNoKHRoaXMuZWRnZSk7XG4gICAgICAgICAgICAgICAgICAgIC8vIGZvdW5kIGEgY2l0eVxuICAgICAgICAgICAgICAgICAgICBpZiAocmVzLmNpdHlSZWFjaGVkICYmICh0aGlzLm9yaWdpbi5pZCAhPSB0aGlzLmxhc3RDaXR5KSApe1xuICAgICAgICAgICAgICAgICAgICAgICAgLy8gY29tcHV0ZSB0aGUgbGVuZ3RoIG9mIHRoZSBwYXRoXG4gICAgICAgICAgICAgICAgICAgICAgICB2YXIgcGF0aExlbmd0aCA9IHRoaXMuZWRnZXMubWFwKGZ1bmN0aW9uKGUpe3JldHVybiBlLmRpc3RhbmNlfSkucmVkdWNlKGZ1bmN0aW9uKGEsYil7cmV0dXJuIGEgKyBifSk7XG4gICAgICAgICAgICAgICAgICAgICAgICB2YXIgZGVsdGFQaGVyb21vbmUgPSAxL3BhdGhMZW5ndGg7XG4gICAgICAgICAgICAgICAgICAgICAgICB0aGlzLmVkZ2VzLmZvckVhY2goZnVuY3Rpb24oZSl7XG4gICAgICAgICAgICAgICAgICAgICAgICAgICAgdmFyIGEgPSBlLnB0MSwgYiA9IGUucHQyLCB3ZWlnaHQgPSAxOyAgXG4gICAgICAgICAgICAgICAgICAgICAgICAgICAgLy8gaW5jcmVhc2VkIGRyb3BwZWQgcGhlcm9tb25zIGZvciB0ZXh0RWRnZXNcbiAgICAgICAgICAgICAgICAgICAgICAgICAgICBpZiAoKGNpdHlTZXQuaW5kZXhPZihhLmlkKSAhPSAtMSkgJiYgY2l0eVNldC5pbmRleE9mKGIuaWQpICE9IC0xICYmIChNYXRoLmFicyhhLmlkIC0gYi5pZCkgPT0gMSkpXG4gICAgICAgICAgICAgICAgICAgICAgICAgICAge1xuICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICB3ZWlnaHQgKj0gMTA7XG4gICAgICAgICAgICAgICAgICAgICAgICAgICAgfVxuICAgICAgICAgICAgICAgICAgICAgICAgICAgIGUucGhlcm9tb24gKz0gKGRlbHRhUGhlcm9tb25lICogd2VpZ2h0KTtcbiAgICAgICAgICAgICAgICAgICAgICAgIH0pO1xuXG4gICAgICAgICAgICAgICAgICAgICAgICB0aGlzLmVkZ2VzID0gW3RoaXMuZWRnZV07XG4gICAgICAgICAgICAgICAgICAgICAgICB0aGlzLmxhc3RDaXR5ID0gdGhpcy5vcmlnaW4uaWQ7XG4gICAgICAgICAgICAgICAgICAgIH1cbiAgICAgICAgICAgICAgICB9XG4gICAgICAgICAgICAgIGJyZWFrO1xuICAgICAgICAgICAgfVxuXG4gICAgICAgIH0sXG5cbiAgICAgICAgc2V0RGlyZWN0aW9uOiBmdW5jdGlvbigpe1xuICAgICAgICAgICAgdmFyIHBvc3NpYmxlRWRnZXMgPSBbXTtcblxuICAgICAgICAgICAgZm9yICh2YXIgaSA9IDA7IGkgPCB0aGlzLm9yaWdpbi5uZXh0cy5sZW5ndGg7IGkrKylcbiAgICAgICAgICAgIHtcbiAgICAgICAgICAgICAgICBwb3NzaWJsZUVkZ2VzW2ldID0gdGhpcy5vcmlnaW4ubmV4dHNbaV07XG4gICAgICAgICAgICB9IFxuXG4gICAgICAgICAgICBwb3NzaWJsZUVkZ2VzLnNwbGljZShwb3NzaWJsZUVkZ2VzLmluZGV4T2YodGhpcy5lZGdlKSwxKTtcblxuICAgICAgICAgICAgLy8gZmxpcCBhIGNvaW4gYW5kIGVpdGhlciB0YWtlIHRoZSBzbWVsbGllc3QgcGF0aCBvciBhIHJhbmRvbSBvbmVcbiAgICAgICAgICAgIGlmIChyYW5kb20oKSA+IDAuNSl7XG4gICAgICAgICAgICAgICAgdmFyIHNtZWxscyA9IHBvc3NpYmxlRWRnZXMubWFwKGZ1bmN0aW9uKGUpe3JldHVybiBlLnBoZXJvbW9ufSk7XG4gICAgICAgICAgICAgICAgdmFyIGluZGV4ID0gc21lbGxzLmluZGV4T2YoTWF0aC5tYXguYXBwbHkoTWF0aCwgc21lbGxzKSk7XG4gICAgICAgICAgICAgICAgdGhpcy5lZGdlID0gcG9zc2libGVFZGdlc1tpbmRleF07XG4gICAgICAgICAgICB9IFxuICAgICAgICAgICAgZWxzZVxuICAgICAgICAgICAgICAgIHRoaXMuZWRnZSA9IHBvc3NpYmxlRWRnZXNbZmxvb3IocmFuZG9tKCkqcG9zc2libGVFZGdlcy5sZW5ndGgpXTtcblxuICAgICAgICAgICAgLy8gc2V0IHRoZSBkZXN0aW5hdGlvbiBwb2ludCwgYmVpbmcgZWRnZS5wdDEgb3IgZWRnZS5wdDJcbiAgICAgICAgICAgIHRoaXMuZGVzdGluYXRpb24gPSAodGhpcy5vcmlnaW4gPT0gdGhpcy5lZGdlLnB0MSkgPyB0aGlzLmVkZ2UucHQyIDogdGhpcy5lZGdlLnB0MTtcblxuICAgICAgICAgICAgdGhpcy5kaXJlY3Rpb24ueCA9IHRoaXMuZGVzdGluYXRpb24ueCAtIHRoaXMub3JpZ2luLng7IFxuICAgICAgICAgICAgdGhpcy5kaXJlY3Rpb24ueSA9IHRoaXMuZGVzdGluYXRpb24ueSAtIHRoaXMub3JpZ2luLnk7XG5cbiAgICAgICAgICAgIHRoaXMuZGlyZWN0aW9uLm5vcm1hbGl6ZSgpO1xuICAgICAgICB9LFxuXG4gICAgICAgIG1vdmU6IGZ1bmN0aW9uKCl7XG4gICAgICAgICAgICB2YXIgZWRnZUNoYW5nZWQ7XG4gICAgICAgICAgICB2YXIgY2l0eVJlYWNoZWQgPSBmYWxzZTtcblxuICAgICAgICAgICAgdGhpcy5kaXJlY3Rpb24ueCA9IHRoaXMuZGVzdGluYXRpb24ueCAtIHRoaXMueDsgXG4gICAgICAgICAgICB0aGlzLmRpcmVjdGlvbi55ID0gdGhpcy5kZXN0aW5hdGlvbi55IC0gdGhpcy55O1xuICAgICAgICAgICAgdGhpcy5kaXJlY3Rpb24ubm9ybWFsaXplKCk7XG5cbiAgICAgICAgICAgIC8vIG9uIGVkZ2VcbiAgICAgICAgICAgIC8vIGlmICgodGhpcy5wcm9nIDwgdGhpcy5lZGdlLmRpc3RhbmNlKSB8fCAoY2FsY3VsYXRlRGlzdGFuY2UodGhpcywgdGhpcy5kZXN0aW5hdGlvbikgPiAwLjAwMSl7XG4gICAgICAgICAgICBpZiAoKGNhbGN1bGF0ZURpc3RhbmNlKHRoaXMsIHRoaXMuZGVzdGluYXRpb24pID4gMC4wMDUpKXtcblxuICAgICAgICAgICAgICAgIC8vIGEgZGVsdGEgbW92ZW1lbnQgd2lsbCBiZSBhcHBsaWVkIGlmIGNvbGxpc2lvbiB3aXRoIG9ic3RhY2xlIGRldGVjdGVkXG4gICAgICAgICAgICAgICAgdmFyIGRlbHRhID0gdGhpcy5hdm9pZE9ic3RhY2xlKCk7XG5cbiAgICAgICAgICAgICAgICB0aGlzLnggKz0gdGhpcy52ZWxvY2l0eSAqIHRoaXMuZGlyZWN0aW9uLnggKyBkZWx0YS54ICogMC4wMDU7XG4gICAgICAgICAgICAgICAgdGhpcy55ICs9IHRoaXMudmVsb2NpdHkgKiB0aGlzLmRpcmVjdGlvbi55ICsgZGVsdGEueSAqIDAuMDA1O1xuXG4gICAgICAgICAgICAgICAgdGhpcy5wcm9nID0gdGhpcy5jYWxjdWxhdGVQcm9ncmVzc2lvbigpO1xuICAgICAgICAgICAgICAgIC8vIHRoaXMucHJvZyA9IGNhbGN1bGF0ZURpc3RhbmNlKHRoaXMsIHRoaXMub3JpZ2luKTtcbiAgICAgICAgICAgICAgICAvL2NvbnNvbGUubG9nKHRoaXMucHJvZyAvIHRoaXMuZWRnZS5kaXN0YW5jZSk7XG4gICAgICAgICAgICAgICAgXG4gICAgICAgICAgICAgICAgZWRnZUNoYW5nZWQgPSBmYWxzZTtcblxuICAgICAgICAgICAgLy8gb24gdmVydGV4XG4gICAgICAgICAgICB9IGVsc2Uge1xuICAgICAgICAgICAgICAgIHRoaXMuc3RlcCA9IDA7XG4gICAgICAgICAgICAgICAgdGhpcy5wcm9nID0gMDtcbiAgICAgICAgICAgICAgICB0aGlzLm9yaWdpbiA9IHRoaXMuZGVzdGluYXRpb247XG4gICAgICAgICAgICAgICAgdGhpcy54ID0gdGhpcy5vcmlnaW4ueDtcbiAgICAgICAgICAgICAgICB0aGlzLnkgPSB0aGlzLm9yaWdpbi55O1xuXG4gICAgICAgICAgICAgICAgdGhpcy5zZXREaXJlY3Rpb24oKTtcblxuICAgICAgICAgICAgICAgIGNpdHlSZWFjaGVkID0gKGNpdHlTZXQuaW5kZXhPZih0aGlzLm9yaWdpbi5pZCkgIT0gLTEpO1xuICAgICAgICAgICAgICAgIGVkZ2VDaGFuZ2VkID0gdHJ1ZTtcbiAgICAgICAgICAgIH1cbiAgICAgICAgICAgIHJldHVybiB7Y2l0eVJlYWNoZWQ6IGNpdHlSZWFjaGVkLCBlZGdlQ2hhbmdlZDogZWRnZUNoYW5nZWR9O1xuICAgICAgICB9LFxuXG4gICAgICAgIGF2b2lkT2JzdGFjbGU6IGZ1bmN0aW9uKCl7XG4gICAgICAgICAgICB2YXIgZGlzdGFuY2UgPSBjYWxjdWxhdGVEaXN0YW5jZSh0aGlzLCBtb3VzZSk7XG4gICAgICAgIFxuICAgICAgICAgICAgaWYgKGRpc3RhbmNlIDw9IFJFUFVMU0lPTilcbiAgICAgICAgICAgIHtcbiAgICAgICAgICAgICAgICAvLyBpZiAoZGlzdGFuY2VFZGdlID4gMC4wMSl7XG4gICAgICAgICAgICAgICAgLy8gICAgIHRoaXMuZGlyZWN0aW9uLnggPSB0aGlzLmRlc3RpbmF0aW9uLnggLSB0aGlzLng7IFxuICAgICAgICAgICAgICAgIC8vICAgICB0aGlzLmRpcmVjdGlvbi55ID0gdGhpcy5kZXN0aW5hdGlvbi55IC0gdGhpcy55O1xuICAgICAgICAgICAgICAgIC8vICAgICB0aGlzLmRpcmVjdGlvbi5ub3JtYWxpemUoKTtcbiAgICAgICAgICAgICAgICAvLyB9XG5cbiAgICAgICAgICAgICAgICByZXR1cm4ge1xuICAgICAgICAgICAgICAgICAgICAvLyBkZWx0YSBtb3ZlbWVudCBpcyBjb21wb3NlZCBvZiBhIHJlcHVsc2lvbiBkZWx0YSBhbmQgYSBjaXJjdWxhciBkZWx0YSBcbiAgICAgICAgICAgICAgICAgICAgeDogKHRoaXMueCAtIG1vdXNlLngpL2Rpc3RhbmNlICsgKHRoaXMueSAtIG1vdXNlLnkpL2Rpc3RhbmNlICogMSxcbiAgICAgICAgICAgICAgICAgICAgeTogKHRoaXMueSAtIG1vdXNlLnkpL2Rpc3RhbmNlIC0gKHRoaXMueCAtIG1vdXNlLngpL2Rpc3RhbmNlICogMVxuICAgICAgICAgICAgICAgICAgICAvLyB4OiAodGhpcy54IC0gbW91c2UueCkvZGlzdGFuY2UgKyAodGhpcy55IC0gbW91c2UueSkvZGlzdGFuY2UgKiAxICsgMC4xICogdGhpcy5kaXJlY3Rpb24ueCxcbiAgICAgICAgICAgICAgICAgICAgLy8geTogKHRoaXMueSAtIG1vdXNlLnkpL2Rpc3RhbmNlIC0gKHRoaXMueCAtIG1vdXNlLngpL2Rpc3RhbmNlICogMSArIDAuMSAqIHRoaXMuZGlyZWN0aW9uLnkgXG4gICAgICAgICAgICAgICAgfTtcbiAgICAgICAgICAgIH1cbiAgICAgICAgICAgIGVsc2VcbiAgICAgICAgICAgICAgICByZXR1cm4ge3g6MCwgeTowfTtcbiAgICAgICAgfSxcblxuICAgICAgICBjYWxjdWxhdGVQcm9ncmVzc2lvbjogZnVuY3Rpb24oKXtcbiAgICAgICAgICAgIHZhciB2ID0gbmV3IFZlY3Rvcih0aGlzLnggLSB0aGlzLm9yaWdpbi54LCB0aGlzLnkgLSB0aGlzLm9yaWdpbi55KTtcbiAgICAgICAgICAgIHZhciBub3JtID0gdi5ub3JtKCk7XG5cbiAgICAgICAgICAgIHZhciB0aGV0YSA9ICh2LnggKiB0aGlzLmVkZ2UuZGlyZWN0aW9uLnggKyB2LnkgKiB0aGlzLmVkZ2UuZGlyZWN0aW9uLnkpIC8gbm9ybTtcbiAgICAgICAgICAgIHZhciBwcm9nID0gbm9ybSAqIE1hdGguYWJzKHRoZXRhKTtcbiAgICAgICAgICAgIC8vIHJldHVybnMgbGVuZ3RoIG9mIHByb2plY3Rpb24gb24gZWRnZVxuICAgICAgICAgICAgcmV0dXJuIHByb2c7XG4gICAgICAgIH1cblxuICAgIH07XG4gICAgcmV0dXJuIEFudDtcbn1cblxuIiwiJ3VzZSBzdHJpY3QnXG5cbnZhciBkdCA9IHJlcXVpcmUoXCJkZWxhdW5heS10cmlhbmd1bGF0ZVwiKTtcblxudmFyIHJhbmdlID0gcmVxdWlyZSgnLi91dGlsaXRpZXMuanMnKS5yYW5nZTtcblxudmFyIHBvaW50cyA9IHJlcXVpcmUoJy4vaW5pdGlhbGl6ZVBvaW50cy5qcycpLnBvaW50cztcbnZhciB0ZXh0TWVzaCA9IHJlcXVpcmUoJy4vaW5pdGlhbGl6ZVBvaW50cy5qcycpLnRleHRNZXNoO1xuLy92YXIgdGV4dFBvaW50c0lkID0gcmVxdWlyZSgnLi9pbml0aWFsaXplUG9pbnRzLmpzJykudGV4dFBvaW50c0lkO1xudmFyIGNpdHlTZXQgPSByZXF1aXJlKCcuL2luaXRpYWxpemVQb2ludHMuanMnKS5jaXR5U2V0O1xudmFyIG5iUmFuZG9tUG9pbnRzID0gcmVxdWlyZSgnLi9pbml0aWFsaXplUG9pbnRzLmpzJykubmJSYW5kb21Qb2ludHM7XG52YXIgZm9yY2VkRWRnZXMgPSByZXF1aXJlKCcuL2luaXRpYWxpemVQb2ludHMuanMnKS5mb3JjZWRFZGdlcztcblxudmFyIEVkZ2UgPSByZXF1aXJlKCcuL2VkZ2UuanMnKTtcblxuLy8gdHJpYW5ndWxhdGVcbnZhciBjZWxscyA9IGR0KHBvaW50cy5tYXAoZnVuY3Rpb24ocCl7XG4gICAgcmV0dXJuIFtwLngsIHAueV1cbn0pKVxuXG52YXIgZWRnZXMgPSBbXTtcbnZhciBwZXJtdXRhdGlvbnMgPSBbWzAsMV0sIFswLDJdLCBbMSwyXV07XG5cbi8vIGZvcmNlIHRoZSBlZGdlcyBvZiB0aGUgdGV4dCB0byBiZSBlZGdlcyBvZiB0aGUgZ3JhcGhcbmlmICh0ZXh0TWVzaCkge1xuICAgIHJhbmdlKDAsIHBvaW50cy5sZW5ndGggLSBuYlJhbmRvbVBvaW50cykuZm9yRWFjaChmdW5jdGlvbihpZCl7XG4gICAgICAgIHZhciBkaXJlY3RMaW5rID0gZm9yY2VkRWRnZXNbaWRdO1xuICAgICAgICB2YXIgdGV4dEVkZ2UgPSBFZGdlLmNyZWF0ZShwb2ludHNbaWRdLCBwb2ludHNbZGlyZWN0TGlua10pO1xuICAgICAgICBlZGdlcy5wdXNoKHRleHRFZGdlKTtcbiAgICAgICAgcG9pbnRzW2lkXS5uZXh0cy5wdXNoKHRleHRFZGdlKTtcbiAgICAgICAgLy9hZGRUb05leHRFZGdlcyhwb2ludHNbcHRdLCB1bmRlZmluZWQsIFt0ZXh0RWRnZV0sIG5leHRFZGdlcyk7XG4gICAgfSlcbn1cblxuXG5jZWxscy5mb3JFYWNoKGZ1bmN0aW9uKGNlbGwpe1xuICAgXG4gICAgZm9yICh2YXIgaSA9IDA7IGkgPCAzOyArK2kpeyAgLy8gZm9yIGVhY2ggcG9pbnQuaWQgbGlzdGVkIGluIGN1cnJlbnQgY2VsbFxuICAgICAgICB2YXIgcHQgPSBwb2ludHNbY2VsbFtpXV07XG4gICAgICAgIC8vY29uc29sZS5sb2coXCItLS0tLS0tLSBpZDpcIiArIHB0KTtcblxuICAgICAgICBmb3IgKHZhciBqID0gMTsgaiA8IDM7ICsrail7IFxuXG4gICAgICAgICAgICB2YXIgcHRqID0gcG9pbnRzW2NlbGxbKCBpICsgaiApICUgM11dOyAvLyBwaWNrIG9uZSBvZiB0aGUgb3RoZXIgMiBwb2ludHMgb2YgdGhlIGNlbGxcbiAgICAgICAgICAgIC8vY29uc29sZS5sb2coXCItLS0gb3RoZXI6IFwiICsgcHRqKTtcbiAgICAgICAgICAgIHZhciBuZXdFZGdlID0gdW5kZWZpbmVkO1xuXG4gICAgICAgICAgICAvLyBpZiBwdCBhbHJlYWR5IGhhcyBuZXh0RWRnZXMgLi4uXG4gICAgICAgICAgICBpZiAocHQubmV4dHMubGVuZ3RoICE9IDApIHtcbiAgICAgICAgICAgICAgICAvL2NvbnNvbGUubG9nKFwiVHJvdXbDqSAxXCIpO1xuICAgICAgICAgICAgICAgIFxuICAgICAgICAgICAgICAgIC8vIC4uLiBnZXQgdGhlIHBvaW50cyBjb3JyZXNwb25kaW5nIC4uLlxuICAgICAgICAgICAgICAgIHZhciB0ZW1wUG9pbnRzID0gcHQubmV4dHMubWFwKGZ1bmN0aW9uKGUpe1xuICAgICAgICAgICAgICAgICAgICByZXR1cm4gW2UucHQxLCBlLnB0Ml07XG4gICAgICAgICAgICAgICAgfSkucmVkdWNlKGZ1bmN0aW9uKGEsIGIpe1xuICAgICAgICAgICAgICAgICAgICAgcmV0dXJuIGEuY29uY2F0KGIpO1xuICAgICAgICAgICAgICAgIH0pO1xuXG4gICAgICAgICAgICAgICAgLy9jb25zb2xlLmxvZyh0ZW1wUG9pbnRzKTtcblxuICAgICAgICAgICAgICAgIC8vIC4uLiBhbmQgY2hlY2sgaWYgcHRqIGFscmVhZHkgaXMgcGFydCBvZiB0aGUgZXhpc3RpbmcgbmV4dEVkZ2VzLiBJZiBub3QsIGFkZCB0aGUgZWRnZS5cbiAgICAgICAgICAgICAgICBpZiAodGVtcFBvaW50cy5pbmRleE9mKHB0aikgPT0gLTEpe1xuICAgICAgICAgICAgICAgICAgICAvL2NvbnNvbGUubG9nKFwiVHJvdXbDqSAyXCIpO1xuICAgICAgICAgICAgICAgICAgICAvL25ld0VkZ2UgPSBjcmVhdGVFZGdlKHBvaW50c1twdF0sIHBvaW50c1twdGpdKTtcbiAgICAgICAgICAgICAgICAgICAgbmV3RWRnZSA9IEVkZ2UuY3JlYXRlKHB0LCBwdGopO1xuICAgICAgICAgICAgICAgICAgICBlZGdlcy5wdXNoKG5ld0VkZ2UpO1xuICAgICAgICAgICAgICAgICAgICBwdC5uZXh0cy5wdXNoKG5ld0VkZ2UpO1xuICAgICAgICAgICAgICAgICAgICAvL2NvbnNvbGUubG9nKGVkZ2VzLmxlbmd0aCk7XG4gICAgICAgICAgICAgICAgfVxuICAgICAgICAgICAgfVxuICAgICAgICAgICAgZWxzZSB7XG4gICAgICAgICAgICAgICAgLy9jb25zb2xlLmxvZyhcIlRyb3V2w6kgM1wiKTtcbiAgICAgICAgICAgICAgICBuZXdFZGdlID0gRWRnZS5jcmVhdGUocHQsIHB0aik7XG4gICAgICAgICAgICAgICAgZWRnZXMucHVzaChuZXdFZGdlKTtcbiAgICAgICAgICAgICAgICBwdC5uZXh0cy5wdXNoKG5ld0VkZ2UpO1xuICAgICAgICAgICAgICAgIC8vY29uc29sZS5sb2coZWRnZXMubGVuZ3RoKTtcbiAgICAgICAgICAgIH1cblxuICAgICAgICAgICAgLy8gYWRkIGFsc28gdGhlIGVkZ2UgdG8gdGhlIGVkZ2UncyBvdGhlciBwb2ludCdzIG5leHRFZGdlc1xuICAgICAgICAgICAgaWYgKG5ld0VkZ2UgIT0gdW5kZWZpbmVkKXtcbiAgICAgICAgICAgICAgICBwdGoubmV4dHMucHVzaChuZXdFZGdlKTtcbiAgICAgICAgICAgICAgICAvL2FkZFRvTmV4dEVkZ2VzKHBvaW50c1twdGpdLCBjZWxsLCB0ZW1wLCBuZXh0RWRnZXMpOyBcbiAgICAgICAgICAgIH0gICAgICAgICBcbiAgICAgICAgfVxuXG4gICAgICAgIC8vIGFkZCB0aGUgdGV4dEVkZ2VzIHRvIG5leHRFZGdlcyBtYXBcbiAgICAgICAgLy9pZiAodGV4dE1lc2ggJiYgdGV4dFBvaW50c0lkLmluZGV4T2YocHQuaWQpICE9IC0xICYmIHB0LmlkIDwgKHRleHRQb2ludHNJZC5sZW5ndGggLSAxKSkgXG4gICAgICAgIC8vIGlmICh0ZXh0TWVzaCAmJiBjaXR5U2V0LmhhcyhwdCkpXG4gICAgICAgIGlmICh0ZXh0TWVzaCAmJiAoY2l0eVNldC5pbmRleE9mKHB0LmlkKSAhPSAtMSkpXG4gICAgICAgIHtcbiAgICAgICAgICAgIC8vIGNvbnNvbGUubG9nKCd2ZXJpZicpO1xuICAgICAgICAgICAgdmFyIHRleHRFZGdlID0gRWRnZS5jcmVhdGUocHQsIHBvaW50c1twdC5pZCArIDFdKTtcbiAgICAgICAgICAgIC8vdmFyIHRlbXAyID0gW107XG4gICAgICAgICAgICBlZGdlcy5wdXNoKHRleHRFZGdlKTtcbiAgICAgICAgICAgIHB0Lm5leHRzLnB1c2godGV4dEVkZ2UpO1xuICAgICAgICAgICAgLy90ZW1wMi5wdXNoKHRleHRFZGdlKTtcbiAgICAgICAgICAgIC8vYWRkVG9OZXh0RWRnZXMocG9pbnRzW3B0KzFdLCB1bmRlZmluZWQsIHRlbXAyLCBuZXh0RWRnZXMpO1xuICAgICAgICB9XG5cbiAgICAgICAgLy9jb25zb2xlLmxvZyhcImxvbmd1ZXVyOiBcIiArIHBvaW50c1twdF0ubmV4dHMubGVuZ3RoKTtcbiAgICAgICAgLy9hZGRUb05leHRFZGdlcyhwb2ludHNbcHRdLCBjZWxsLCBuZXh0cywgbmV4dEVkZ2VzKTtcbiAgICB9XG4gICAgLy9jb25zb2xlLmxvZyhlZGdlcy5sZW5ndGgpO1xufSlcblxubW9kdWxlLmV4cG9ydHMgPSBlZGdlczsiLCIndXNlIHN0cmljdCc7XG5cbnZhciBzcXJ0ID0gTWF0aC5zcXJ0O1xudmFyIHBvdyA9IE1hdGgucG93O1xudmFyIGFicyA9IE1hdGguYWJzO1xudmFyIGF0YW4gPSBNYXRoLmF0YW47XG5cbnZhciBWZWN0b3IgPSByZXF1aXJlKCcuL3ZlY3Rvci5qcycpO1xuXG5cbmZ1bmN0aW9uIEVkZ2UocHRBLCBwdEIpIHtcbiAgICB2YXIgZGlzdGFuY2UgPSBzcXJ0KCBwb3cocHRBLnggLSBwdEIueCwgMikgKyBwb3cocHRBLnkgLSBwdEIueSwgMikgKTtcblxuICAgIC8vIGZpbmQgbGluZSBlcXVhdGlvbiBheCArIGJ5ICsgYyA9IDBcbiAgICB2YXIgYSA9IDE7XG4gICAgdmFyIGIgPSAtIChwdEIueCAtIHB0QS54KSAvIChwdEIueSAtIHB0QS55KTtcblxuICAgIC8vIG9yaWVudGF0ZSB2ZWN0b3IgKGEsYilcbiAgICBpZiAoYiA8IDApe1xuICAgICAgICBiID0gLWI7XG4gICAgICAgIGEgPSAtYTtcbiAgICB9XG5cbiAgICAvLyBub3JtYWxpemUgdmVjdG9yIChhLGIpXG4gICAgdmFyIG4gPSBuZXcgVmVjdG9yKGEsIGIpO1xuICAgIG4ubm9ybWFsaXplKCk7XG5cbiAgICB2YXIgYyA9IC0gKGEgKiBwdEEueCArIGIgKiBwdEEueSk7XG5cbiAgICAvLyAvLyBjYWxjdWxhdGUgdmVjdG9yIGRpcmVjdG9yXG4gICAgdmFyIHYgPSBuZXcgVmVjdG9yKHB0Qi54IC0gcHRBLngsIHB0Qi55IC0gcHRBLnkpO1xuICAgIFxuICAgIHYubm9ybWFsaXplKCk7XG5cblxuICAgIHRoaXMuaWQgPSB1bmRlZmluZWQ7XG4gICAgdGhpcy5wdDEgPSBwdEE7XG4gICAgdGhpcy5wdDIgPSBwdEI7XG4gICAgdGhpcy5kaXJlY3Rpb24gPSB2O1xuICAgIHRoaXMub3J0aERpcmVjdGlvbiA9IG47IFxuICAgIHRoaXMuZGlzdGFuY2UgPSBkaXN0YW5jZTtcbiAgICB0aGlzLnBoZXJvbW9uID0gMS9kaXN0YW5jZTtcbiAgICB0aGlzLmxpbmUgPSB7XG4gICAgICAgIGE6IGEsXG4gICAgICAgIGI6IGIsXG4gICAgICAgIGM6IGMsXG4gICAgfVxufVxuXG5cbi8vIHN0YXRpYyBtZXRob2RzXG5FZGdlLmNyZWF0ZSA9IGZ1bmN0aW9uKHB0QSwgcHRCKSB7XG4gICAgdmFyIGVkZ2UgPSBuZXcgRWRnZShwdEEsIHB0Qik7XG4gICAgcmV0dXJuIGVkZ2U7XG59XG5cblxuLy8gbWV0aG9kc1xuRWRnZS5wcm90b3R5cGUgPSB7XG5cbiAgICBnZXRPdGhlclBvaW50OiBmdW5jdGlvbihwb2ludCkge1xuICAgICAgICBpZiAocG9pbnQgPT0gdGhpcy5wdDEpXG4gICAgICAgICAgICByZXR1cm4gdGhpcy5wdDI7XG4gICAgICAgIGVsc2UgaWYgKHBvaW50ID09IHRoaXMucHQyKVxuICAgICAgICAgICAgcmV0dXJuIHRoaXMucHQxO1xuICAgICAgICBlbHNlXG4gICAgICAgICAgICBjb25zb2xlLmxvZyhcIkVycm9yXCIpO1xuICAgIH0sXG5cbiAgICBjYWxjdWxhdGVEaXN0YW5jZTogZnVuY3Rpb24oeCwgeSkge1xuICAgICAgICB2YXIgYSA9IHRoaXMubGluZS5hLFxuICAgICAgICAgICAgYiA9IHRoaXMubGluZS5iLFxuICAgICAgICAgICAgYyA9IHRoaXMubGluZS5jO1xuICAgICAgICByZXR1cm4gYWJzKGEgKiB4ICsgYiAqIHkgKyBjKSAvIE1hdGguc3FydChNYXRoLnBvdyhhLDIpICsgTWF0aC5wb3coYiwyKSk7XG4gICAgfSxcblxuICAgIC8qY2FsY3VsYXRlRGlzdGFuY2U6IGZ1bmN0aW9uKHBvaW50KSB7XG4gICAgICAgIHZhciBhID0gdGhpcy5saW5lLmEsXG4gICAgICAgICAgICBiID0gdGhpcy5saW5lLmIsXG4gICAgICAgICAgICBjID0gdGhpcy5saW5lLmM7XG4gICAgICAgIHJldHVybiBhYnMoYSAqIHBvaW50LnggKyBiICogcG9pbnQueSArIGMpIC8gTWF0aC5zcXJ0KE1hdGgucG93KGEsMikgKyBNYXRoLnBvdyhiLDIpKTtcbiAgICB9Ki9cbn1cbm1vZHVsZS5leHBvcnRzID0gRWRnZTsiLCIndXNlIHN0cmljdCdcblxudmFyIGFudEZ1bmN0aW9uID0gcmVxdWlyZSgnLi9hbnQuanMnKTtcblxudmFyIG5iQW50cyA9IDMwMDA7XG5cbm1vZHVsZS5leHBvcnRzID0gZnVuY3Rpb24gKGNvbnRhaW5lcikge1xuXG5cdHZhciBBbnQgPSBhbnRGdW5jdGlvbihjb250YWluZXIpO1xuXG5cdHZhciBwb3B1bGF0aW9uID0gbmV3IEFycmF5KG5iQW50cyk7XG5cdHZhciBwb3NzaWJsZVN0YXJ0UG9pbnRzSWQgPSByZXF1aXJlKCcuL2luaXRpYWxpemVQb2ludHMuanMnKS5wb3NzaWJsZVN0YXJ0UG9pbnRzSWQ7XG5cblx0Zm9yICh2YXIgaSA9IDA7IGkgPCBuYkFudHM7IGkrKykge1xuXHQgICAgdmFyIG5ld0FudCA9IG5ldyBBbnQoQW50LmdlbmVyYXRlUmFuZFN0YXJ0UG9pbnQoKSk7XG5cdCAgICBuZXdBbnQuc2V0RGlyZWN0aW9uKCk7XG5cdCAgICBwb3B1bGF0aW9uW2ldID0gbmV3QW50O1xuXHR9XG5cblx0cmV0dXJuIHBvcHVsYXRpb247XG5cbn0iLCIndXNlIHN0cmljdCdcblxudmFyIHBhcnNlID0gcmVxdWlyZSgncGFyc2Utc3ZnLXBhdGgnKTtcblxudmFyIHJhbmdlID0gcmVxdWlyZSgnLi91dGlsaXRpZXMuanMnKS5yYW5nZTtcblxudmFyIFBvaW50ID0gcmVxdWlyZSgnLi9wb2ludC5qcycpO1xuXG52YXIgcmFuZG9tID0gTWF0aC5yYW5kb207XG5cbnZhciBuYlJhbmRvbVBvaW50cyA9IDMwMDtcbnZhciBuYlN0YXJ0UG9pbnRzID0gMjA7XG5cbnZhciBuYkNpdHkgPSAyO1xuXG52YXIgdGV4dE1lc2ggPSB0cnVlO1xuXG4vLyBGcmFtZSBkZWZpbml0aW9uXG52YXIgeEluaXQgPSAwLCB5SW5pdCA9IDA7XG52YXIgdyA9IDEsXG4gICAgaCA9IDE7XG5cbnZhciBBY2hhciA9IFwiYyA0LjYwMTEsLTExLjcxMDQ3IDkuMjA4MzUsLTIzLjQyMDA2IDEzLjgxOTksLTM1LjEyODk4IDQuNjExNTYsLTExLjcwODkyIDkuMjI3NDEsLTIzLjQxNzE4IDEzLjg0NTczLC0zNS4xMjUgNC42MTgzMSwtMTEuNzA3ODIgOS4yMzkwOCwtMjMuNDE1MTkgMTMuODYwNDYsLTM1LjEyMjMzIDQuNjIxMzgsLTExLjcwNzE0IDkuMjQzMzYsLTIzLjQxNDA2IDEzLjg2NDExLC0zNS4xMjA5NyA0LjYyMDc0LC0xMS43MDY5MSA5LjI0MDI1LC0yMy40MTM4MSAxMy44NTY2NywtMzUuMTIwOTIgNC42MTY0MSwtMTEuNzA3MTIgOS4yMjk3NCwtMjMuNDE0NDQgMTMuODM4MTQsLTM1LjEyMjE4IDQuNjA4MzksLTExLjcwNzc1IDkuMjExODYsLTIzLjQxNTkyIDEzLjgwODU0LC0zNS4xMjQ3NCA0LjU5NjY4LC0xMS43MDg4MSA5LjE4NjU4LC0yMy40MTgyNzcgMTMuNzY3ODUsLTM1LjEyODYwMyAxMi45OTkyMywtMy4zNTc4NTUgMjQuMzAwNjksLTYuODIxMTQ0IDMzLjg2ODkzLC00LjQ2NDExIDkuNTY4MjUsMi4zNTcwMzUgMTcuNDAzMjgsMTAuNTM0MzkzIDIzLjQ2OTY3LDMwLjQ1NzgzMyA0LjE3NTk4LDEwLjYyMTE0IDguMzYwMzEsMjEuMjM4ODIgMTIuNTQ5NDIsMzEuODU0NTMgNC4xODkxLDEwLjYxNTcxIDguMzgyOTgsMjEuMjI5NDMgMTIuNTc4MDcsMzEuODQyNjYgNC4xOTUwOSwxMC42MTMyMyA4LjM5MTQsMjEuMjI1OTUgMTIuNTg1MzUsMzEuODM5NjUgNC4xOTM5NSwxMC42MTM2OSA4LjM4NTU1LDIxLjIyODM2IDEyLjU3MTIzLDMxLjg0NTQ5IDQuMTg1NjgsMTAuNjE3MTIgOC4zNjU0NCwyMS4yMzY3IDEyLjUzNTcyLDMxLjg2MDIxIDQuMTcwMjgsMTAuNjIzNSA4LjMzMTA4LDIxLjI1MDkzIDEyLjQ3ODgzLDMxLjg4Mzc3IDQuMTQ3NzUsMTAuNjMyODMgOC4yODI0NSwyMS4yNzEwNyAxMi40MDA1MywzMS45MTYyIDQuMTE4MDgsMTAuNjQ1MTIgOC4yMTk1NSwyMS4yOTcxMiAxMi4zMDA4NSwzMS45NTc0OSAtNS45MDg5Niw2Ljk1NTYxIC0yNC42MTYxNywxLjExMjk4IC0zNS41OTM3MiwzIC0xNi43NTI0OCwyLjg0ODU5IC0yNi45NjQyMSwtMC40MTQxNiAtMjguNDA2MjgsLTE5IC0zLjE3NzI2LC04LjQyMTU3IC02LjM1NjA2LC0xNi44NzkyNCAtOS40Nzk5NywtMjUuMzQ4NTQgLTMuMTIzOTEsLTguNDY5MjkgLTYuMTkyOTMsLTE2Ljk1MDIgLTkuMTUwNjIsLTI1LjQxODI2IC0xNi40NjcyMywtMC44MDA1MyAtMzUuMTc3NjgsLTIuNzQyODEgLTUzLjMzNzI3LC0zLjExNDM5IC0xOC4xNTk2LC0wLjM3MTU4IC0zNS43NjgzNCwwLjgyNzUzIC01MC4wMzIxNCw2LjMwOTc2IC00LjE1Njc5LDEwLjkzMTI0IC04LjEzODA2LDIxLjkyNjkgLTEyLjE1Nzg1LDMyLjkwNTQ2IC00LjAxOTgsMTAuOTc4NTUgLTguMDc4MTMsMjEuOTQwMDEgLTEyLjM4OTAzLDMyLjgwMjgyIC05LjUyNzU4LDAuODI3NDcgLTE5LjEwNDU3LDAuOTY0MjMgLTI4LjY5MjgyLDAuOTMzNjMgLTkuNTg4MjQsLTAuMDMwNiAtMTkuMTg3NzMsLTAuMjI4NTUgLTI4Ljc2MDMsLTAuMDcwNSBsIDAsLTIuMTk3OTEgeiBtIDE3NCwtMTA5IGMgLTMuMTc0NjIsLTkuNTExMzYgLTYuNDQ4MzUsLTE4Ljk5MTc0IC05Ljc1NDMsLTI4LjQ2MjI0IC02LjgxNTU5LC0xOS41MTQxMiAtNi4yNDIwMiwtMjUuMzc5NDYgLTEyLjMxNzQsLTQzLjYwNzg4IC0zLjE2NzIyLC05LjUxNTQzIC0xMy42MDI1NywtMzIuMzI4NjYgLTE2LjQ5OTczLC00MS45Mjk4OCAtNi4wODk4OSwxMi44ODU3NiAtMTEuMTMyLDI2LjU5MjcyIC0xNS45MzQxNSw0MC40NzQ3NiAtNC44MDIxNSwxMy44ODIwMyAtOS4zNjQzNSwyNy45MzkxNSAtMTQuNDk0NDIsNDEuNTI1MjQgLTIuMzk0MjUsMTIuMDU4MDEgLTIyLjg4MTUsNDAuMTExMTYgMi4xODc0NSwzNCAxMS4wNTI1NiwtMC40NDQ4NiAyMi41MjkzOSwwLjExMDYxIDMzLjg1NjI0LDAuMjQ5NTUgMTEuMzI2ODQsMC4xMzg5NSAyMi41MDM2OSwtMC4xMzg2MyAzMi45NTYzMSwtMi4yNDk1NSB6IFwiO1xudmFyIE5jaGFyID0gXCJjIDAsMCAwLC0yMy44MzMzNCAwLC0zNS43NSAwLC0xMS45MTY2NyAwLC0yMy44MzMzNCAwLC0zNS43NSAwLC0xMS45MTY2NiAwLC0yMy44MzMzMyAwLC0zNS43NSAwLC0xMS45MTY2NjcgMCwtMjMuODMzMzM1IDAsLTM1Ljc1MDAwMyAxNC42NDcyOSwxLjA1MzEzIDMxLjAzMTkzLC0yLjA4ODA4IDQ0LjYwNjc5LDEuNTU0MTUgNi4xODA3Niw3LjYxMDk5NiAxMi4zNTQzOCwxNS4yMzEzNDIgMTguNTE5NzMsMjIuODYxMDEzIDYuMTY1MzUsNy42Mjk2NyAxMi4zMjI0NCwxNS4yNjg2NiAxOC40NzAxNCwyMi45MTY5NiA2LjE0NzcsNy42NDgyOSAxMi4yODYwMiwxNS4zMDU4OCAxOC40MTM4MywyMi45NzI3NCA2LjEyNzgxLDcuNjY2ODYgMTIuMjQ1MTIsMTUuMzQzIDE4LjM1MDgxLDIzLjAyODM4IDYuMTA1NjgsNy42ODUzOCAxMi4xOTk3NSwxNS4zODAwMSAxOC4yODEwNywyMy4wODM4NiA2LjA4MTMyLDcuNzAzODQgMTIuMTQ5OSwxNS40MTY5MSAxOC4yMDQ2MiwyMy4xMzkxOCA2LjA1NDcyLDcuNzIyMjYgMTIuMDk1NTcsMTUuNDUzNzIgMTguMTIxNDUsMjMuMTk0MzUgNi4wMjU4Nyw3Ljc0MDYzIDEyLjAzNjc2LDE1LjQ5MDQzIDE4LjAzMTU2LDIzLjI0OTM3IDAuNDU1OCwtMTUuNDkxNCAwLjcyNjU1LC0zMC45ODcyNSAwLjg4MTE5LC00Ni40ODU4MiAwLjE1NDY0LC0xNS40OTg1OCAwLjE5MzE5LC0zMC45OTk4OCAwLjE4NDU4LC00Ni41MDIxOCAtMC4wMDksLTE1LjUwMjMgLTAuMDY0MywtMzEuMDA1NiAtMC4wOTgzLC00Ni41MDgxNyAtMC4wMzQsLTE1LjUwMjU4IC0wLjA0NjEsLTMxLjAwNDQzMiAwLjAzMjUsLTQ2LjUwMzgzMyA5LjMzMzMzLDAgMTguNjY2NjcsMCAyOCwwIDkuMzMzMzMsMCAxOC42NjY2NiwwIDI4LDAgMCwxMS45MTY2NjcgMCwyMy44MzMzMzIgMCwzNS43NTAwMDMgMCwxMS45MTY2NiAwLDIzLjgzMzMzIDAsMzUuNzUgMCwxMS45MTY2NiAwLDIzLjgzMzMzIDAsMzUuNzUgMCwxMS45MTY2NyAwLDIzLjgzMzMzIDAsMzUuNzUgMCwxMS45MTY2NiAwLDIzLjgzMzMzIDAsMzUuNzUgMCwxMS45MTY2NyAwLDIzLjgzMzM0IDAsMzUuNzUgMCwxMS45MTY2NiAwLDIzLjgzMzMzIDAsMzUuNzUgMCwxMS45MTY2NyAwLDIzLjgzMzM0IDAsMzUuNzUgLTE0LjY0NDEsLTEuMDUzNDggLTMxLjAyNiwyLjA4ODQ3IC00NC41OTc2NCwtMS41NTQxNiAtNi4yNzI0NywtNy4zODAzNyAtMTIuNDg0MzYsLTE0LjgxNjkyIC0xOC42NTEyMywtMjIuMjk1MDcgLTYuMTY2ODgsLTcuNDc4MTUgLTEyLjI4ODczLC0xNC45OTc5IC0xOC4zODExMSwtMjIuNTQ0NjUgLTYuMDkyMzgsLTcuNTQ2NzUgLTEyLjE1NTI4LC0xNS4xMjA1MSAtMTguMjA0MjUsLTIyLjcwNjY3IC02LjA0ODk2LC03LjU4NjE3IC0xMi4wODM5OSwtMTUuMTg0NzYgLTE4LjEyMDYzLC0yMi43ODExNiAtNi4wMzY2NCwtNy41OTY0IC0xMi4wNzQ4OSwtMTUuMTkwNjIgLTE4LjEzMDMsLTIyLjc2ODA3IC02LjA1NTQxLC03LjU3NzQ1IC0xMi4xMjc5NywtMTUuMTM4MTMgLTE4LjIzMzIzLC0yMi42Njc0NCAtNi4xMDUyNiwtNy41MjkzMiAtMTIuMjQzMjIsLTE1LjAyNzI3IC0xOC40Mjk0MiwtMjIuNDc5MjYgLTYuMTg2MTksLTcuNDUxOTkgLTEyLjQyMDYzLC0xNC44NTgwMyAtMTguNzE4ODYsLTIyLjIwMzUyIC0wLjIyNzkxLDE1LjE2NDk2IC0wLjM2Njc0LDMwLjMzMDggLTAuNDQ4OSw0NS40OTcxOCAtMC4wODIyLDE1LjE2NjM3IC0wLjEwNzY2LDMwLjMzMzI5IC0wLjEwODg5LDQ1LjUwMDQgLTAuMDAxLDE1LjE2NzEyIDAuMDIxOCwzMC4zMzQ0MyAwLjAzNjcsNDUuNTAxNjIgMC4wMTQ5LDE1LjE2NzE4IDAuMDIxNiwzMC4zMzQyMiAtMC4wMTIyLDQ1LjUwMDggLTkuMzMzMzMsMCAtMTguNjY2NjcsMCAtMjgsMCAtOS4zMzMzMywwIC0xOC42NjY2NiwwIC0yOCwwIDAsLTExLjkxNjY3IDAsLTIzLjgzMzM0IDAsLTM1Ljc1IDAsLTExLjkxNjY3IDAsLTIzLjgzMzM0IDAsLTM1Ljc1IDAsLTExLjkxNjY2IDAsLTIzLjgzMzMzIDAsLTM1Ljc1IDAsLTExLjkxNjY3IDAsLTM1Ljc1IC0xMGUtNiwtMzUuNzUgeiBcIjtcbnZhciBUY2hhciA9IFwiYyAwLC05LjkxNjY3IDAsLTE5LjgzMzM0IDAsLTI5Ljc1IDAsLTkuOTE2NjcgMCwtMTkuODMzMzQgMCwtMjkuNzUgMCwtOS45MTY2NiAwLC0xOS44MzMzMyAwLC0yOS43NSAwLC05LjkxNjY2IC0xLjQ3Mzk0LC0xOS44MzMzMyAtMS40NzM5NCwtMjkuNzUgMCwtMTUuMzMzMzQgLTI5LjE5MjczLDAgLTQ0LjUyNjA2LDAgLTE1LjMzMzMzLDAgLTMwLjY2NjY3LDAgLTQ2LDAgMCwtOCAwLC0xNiAwLC0yNC4wMDAwMDIgMCwtOCAwLC0xNi4wMDAwMDEgMCwtMjQuMDAwMDAxIDkuOTE2NjYsMCAxOS44MzMzMywwIDI5Ljc1LDAgOS45MTY2NywwIDE5LjgzMzMzLDAgMjkuNzUsMCA5LjkxNjY2LDAgMTkuODMzMzMsMCAyOS43NSwwIDkuOTE2NjcsMCAxOS44MzMzMywwIDI5Ljc1LDAgOS45MTY2NiwwIDE5LjgzMzMzLDAgMjkuNzUsMCA5LjkxNjY3LDAgMTkuODMzMzMsMCAyOS43NSwwIDkuOTE2NjYsMCAxOS44MzMzMywwIDI5Ljc1LDAgOS45MTY2NywwIDE5LjgzMzMzLDAgMjkuNzUsMCAwLDggMCwxNi4wMDAwMDEgMCwyNC4wMDAwMDEgMCw4LjAwMDAwMiAwLDE2LjAwMDAwMiAwLDI0LjAwMDAwMiAtMTUsMCAtMzAsMCAtNDUsMCAtMTUsMCAtNDQuODg4MDksLTEzLjUyNTY1IC00NSwxLjQ3Mzk0IC0wLjA3MzUsOS44NTU4OCAtMC4xMDYxMywxOC4yMzk5IC0wLjExMjQ5LDI4LjA5OTIzIC0wLjAwNiw5Ljg1OTMyIDAuMDEzNSwxOS43MjAwMSAwLjA0NSwyOS41ODEzMyAwLjAzMTUsOS44NjEzMiAwLjA3NDUsMTkuNzIzMjkgMC4xMTQzLDI5LjU4NTE3IDAuMDM5OSw5Ljg2MTg4IDAuMDc2NSwxOS43MjM2NyAwLjA5NTQsMjkuNTg0NjcgMC4wMTg5LDkuODYwOTkgMC4wMTk5LDE5LjcyMTE4IC0wLjAxMTcsMjkuNTc5ODQgLTAuMDMxNSw5Ljg1ODY2IC0wLjA5NTYsMTkuNzE1OCAtMC4yMDY4OCwyOS41NzA2OSAtMC4xMTEzMSw5Ljg1NDg4IC0wLjI2OTg1LDE5LjcwNzUyIC0wLjQ5MDMzLDI5LjU1NzE5IC0wLjIyMDQ3LDkuODQ5NjcgLTAuNTAyODksMTkuNjk2MzYgLTAuODYxOTMsMjkuNTM5MzcgLTcuMzk2MjQsLTAuOTk5NjQgLTE4Ljk1NjU4LDAuOTY3MjYgLTI5LjcwOTEyLDEuODI5NzEgLTEwLjc1MjUzLDAuODYyNDUgLTI0LjQ4NTU2LDIuMDI2MDkgLTI0Ljg2MjMxLC00Ljc5Njk2IC0wLjY2MjIzLC0xMS45OTMxNCAtMC42NjIyMywtMjEuNTQzNDkgLTAuNDk2NjcsLTMwLjQ4MzE0IDAuMTY1NTYsLTguOTM5NjUgMC40OTY2NywtMTcuMjY4NiAwLjQ5NjY3LC0yNi44MTg5NSAwLC05LjU1MDM1IDAsLTE5LjEwMDcgMCwtMjguNjUxMDUgMCwtOS41NTAzNSAwLC0xOS4xMDA2OSAwLC0yOC42NTEwNCB6XCI7XG52YXIgU2NoYXIgPSBcImMgLTE2LjY1ODcsLTIuMzQzODcgLTMzLjM1OTk1LC02LjIzNTE1IC00OS4yNTEzNywtMTIuMDg4MzEgLTE1Ljg5MTQyLC01Ljg1MzE2IC0zMC45NzMwNSwtMTMuNjY4MiAtNDQuMzkyNiwtMjMuODU5NTcgMy43MzcxNSwtOC4wMjYzOSA3Ljk2MDY1LC0xNS43OTM3MSAxMi4xOTEwMSwtMjMuNTUyMjggNC4yMzAzNywtNy43NTg1OCA4LjQ2NzU5LC0xNS41MDg0MiAxMi4yMzIxOSwtMjMuNDk5ODQgOC44MjkwMiw2LjY5MTYyIDE4LjQ4MjM1LDEyLjc1NjYzIDI4LjY3NDEyLDE3Ljg2MzQ0IDEwLjE5MTc4LDUuMTA2ODEgMjAuOTIyMDUsOS4yNTU0MiAzMS45MDQ4NSwxMi4xMTQyMSAxMC45ODI5LDIuODU4OCAyMi4yMTgzLDQuNDI3OCAzMy40MjA2LDQuMzc1MzkgMTEuMjAyMiwtMC4wNTI0IDIyLjM3MTMsLTEuNzI2MjIgMzMuMjIxMiwtNS4zNTMwNCAxNC4yODA0LC02LjkyMTY5IDE3LjMwMzMsLTE5LjY0NDM2IDEzLjQ5NDYsLTMxLjE1MDI4IC0zLjgwODcsLTExLjUwNTkxIC0xNC40NDg5LC0yMS43OTUwNyAtMjcuNDk0NiwtMjMuODQ5NzIgLTEwLjE2MiwtNC4wNzY4NiAtMjEuMTA1MiwtNy4xNDg2MyAtMzIuMTk3NSwtMTAuMTEwOTcgLTExLjA5MjQsLTIuOTYyMzUgLTIyLjMzNCwtNS44MTUyNiAtMzMuMDkzMSwtOS40NTQzOSAtMTAuNzU5MSwtMy42MzkxNCAtMjEuMDM1NjIsLTguMDY0NDkgLTMwLjE5Nzc5LC0xNC4xNzE3MSAtOS4xNjIxNiwtNi4xMDcyMyAtMTcuMjA5OTYsLTEzLjg5NjMyIC0yMy41MTE1OCwtMjQuMjYyOTMgLTUuMTc5NCwtMTAuOTcwNTggLTcuNDg1NDQsLTIyLjc3NTgzIC03LjMwNTEyLC0zNC41Mjk2NiAwLjE4MDMxLC0xMS43NTM4MiAyLjg0Njk4LC0yMy40NTYyIDcuNjEyOTksLTM0LjIyMTAzIDQuNzY2MDEsLTEwLjc2NDgzIDExLjYzMTM3LC0yMC41OTIxIDIwLjIwOTA1LC0yOC41OTU2OTIgOC41Nzc2OSwtOC4wMDM1OTIgMTguODY3NzEsLTE0LjE4MzUwNiAzMC40ODMwNSwtMTcuNjUzNjIxIDEyLjE4MjEsLTQuMDI0MTQ1IDI0Ljg3NzcsLTYuMzU3MDA5IDM3LjY5ODgsLTcuMTIxOTAxIDEyLjgyMTEsLTAuNzY0ODkyIDI1Ljc2NzcsMC4wMzgxOSAzOC40NTE3LDIuMjg1OTMzIDEyLjY4NDEsMi4yNDc3NDQgMjUuMTA1NSw1Ljk0MDE1MyAzNi44NzY1LDEwLjk1MzkxOCAxMS43NzA5LDUuMDEzNzY0IDIyLjg5MTIsMTEuMzQ4ODg1IDMyLjk3MywxOC44ODIwNTMgLTQuNDc0NSw3LjAwMzY4IC04LjI1ODEsMTQuNzA3MDggLTEyLjE4ODYsMjIuMjM0OTkgLTMuOTMwNiw3LjUyNzkgLTguMDA4LDE0Ljg4MDMxIC0xMy4wNzAxLDIxLjE4MTk4IC0xNS4wNzI4LC0xMS4xNTE1OCAtMzMuNTAwOSwtMjAuNTQ0OTEgLTUyLjYxMTQsLTI0LjkxNzI2IC0xOS4xMTA1LC00LjM3MjM1IC0zOC45MDM0LC0zLjcyMzcyIC01Ni43MDU4LDUuMjA4NjEgLTEwLjAzMjEsNy4wMjc4OSAtMTMuMzM0OCwxOC44NDY5NSAtMTEuMTU2MSwyOS41MzU5NiAyLjE3ODcsMTAuNjg5MDIgOS44Mzg3LDIwLjI0Nzk5IDIxLjczMiwyMi43NTU3MiAxMC4xMzM1LDQuNDMzMDcgMjEuMDI0MSw3LjY5NDA3IDMyLjExMDEsMTAuNjkwMjMgMTEuMDg2MSwyLjk5NjE2IDIyLjM2NzYsNS43Mjc0NyAzMy4yODMsOS4xMDExNyAxMC45MTU1LDMuMzczNyAyMS40NjQ4LDcuMzg5NzggMzEuMDg2NSwxMi45NTU0NyA5LjYyMTcsNS41NjU2OSAxOC4zMTU3LDEyLjY4MDk5IDI1LjUyMDQsMjIuMjUzMTMgNi42MzAxLDkuNzE2MzUgMTAuNTU4NywyMC44ODAyMyAxMi4wMDY0LDMyLjQxMjMxIDEuNDQ3NiwxMS41MzIwNyAwLjQxNDIsMjMuNDMyMzQgLTIuODc5NywzNC42MjE0NiAtMy4yOTM5LDExLjE4OTEyIC04Ljg0ODQsMjEuNjY3MDkgLTE2LjQ0MywzMC4zNTQ1NyAtNy41OTQ2LDguNjg3NDkgLTE3LjIyOTMsMTUuNTg0NDkgLTI4LjY4MzcsMTkuNjExNjYgLTEzLjAwOTEsNS45MjkyOCAtMjcuMDE5Nyw4LjYzODY5IC00MS4yNzI4LDkuNjM2MDggLTE0LjI1MywwLjk5NzM4IC0yOC43NDg0LDAuMjgyNzQgLTQyLjcyNzIsLTAuNjM2MDggeiBcIjtcblxudmFyIHN2Z1N0cmluZyA9IFwibSAwLDAgXCIgKyBBY2hhciArXCJtIDEyNiwtMzEgXCIgKyBOY2hhciArIFwibSAzNzYsMjQgXCIrIFRjaGFyICtcIm0gMjUwLDEyMCBcIisgU2NoYXI7XG5cbiBcbmZ1bmN0aW9uIHN2Z1RvUG9pbnRzKHN2Z1N0cmluZykge1xuICAgIHZhciBwb2ludHMgPSBbXTtcbiAgICAvL3ZhciBlZGdlcyA9IFtdO1xuXG4gICAgLy8gdmFyIGVkZ2VzID0gbmV3IE1hcCgpOyAvLyBwb2ludElkIC0+IHBvaW50SWQgb24gYm9yZGVyXG4gICAgdmFyIGVkZ2VzID0gT2JqZWN0LmNyZWF0ZShudWxsKTtcblxuICAgIHZhciBiZWdpbmluZ1BhdGg7XG5cbiAgICB2YXIgWCA9IDA7XG4gICAgdmFyIFkgPSAwO1xuICAgIHZhciBuYlBvaW50cyA9IDA7XG4gICAgdmFyIHByZXZQb2ludDtcblxuICAgIHZhciBjb21tYW5kcyA9IHBhcnNlKHN2Z1N0cmluZylcbiAgICBmb3IgKHZhciBpPTA7IGk8Y29tbWFuZHMubGVuZ3RoOyBpKyspe1xuICAgICAgICB2YXIgY29tbWFuZCA9IGNvbW1hbmRzW2ldO1xuICAgICAgICBzd2l0Y2ggKGNvbW1hbmRbMF0pIHtcbiAgICAgICAgICAgIGNhc2UgXCJtXCI6XG4gICAgICAgICAgICAgICAgWCArPSBjb21tYW5kWzFdO1xuICAgICAgICAgICAgICAgIFkgKz0gY29tbWFuZFsyXTtcbiAgICAgICAgICAgICAgICBwcmV2UG9pbnQgPSB1bmRlZmluZWQ7XG4gICAgICAgICAgICAgICAgYmVnaW5pbmdQYXRoID0gbmJQb2ludHM7XG4gICAgICAgICAgICAgICAgYnJlYWs7XG4gICAgICAgICAgICBjYXNlIFwiTVwiOlxuICAgICAgICAgICAgICAgIFggPSBjb21tYW5kWzFdO1xuICAgICAgICAgICAgICAgIFkgPSBjb21tYW5kWzJdO1xuICAgICAgICAgICAgICAgIHByZXZQb2ludCA9IHVuZGVmaW5lZDtcbiAgICAgICAgICAgICAgICBiZWdpbmluZ1BhdGggPSBuYlBvaW50cztcbiAgICAgICAgICAgICAgICBicmVhazsgIFxuICAgICAgICAgICAgY2FzZSBcImNcIjpcbiAgICAgICAgICAgICAgICBYICs9IGNvbW1hbmRbNV07XG4gICAgICAgICAgICAgICAgWSArPSBjb21tYW5kWzZdO1xuICAgICAgICAgICAgICAgIHBvaW50cy5wdXNoKHtpZDpuYlBvaW50cywgeDpYLCB5Oll9KTtcbiAgICAgICAgICAgICAgICAvLyBuYlBvaW50cysrO1xuICAgICAgICAgICAgICAgIC8vIGlmIChwcmV2UG9pbnQpIHtcbiAgICAgICAgICAgICAgICAvLyAgICAgZWRnZXMucHVzaChbcHJldlBvaW50LCBuYlBvaW50c10pO1xuICAgICAgICAgICAgICAgIC8vICAgICBwcmV2UG9pbnQgPSBuYlBvaW50cztcbiAgICAgICAgICAgICAgICAvLyB9XG4gICAgICAgICAgICAgICAgaWYgKHByZXZQb2ludCAhPSB1bmRlZmluZWQpIHtcbiAgICAgICAgICAgICAgICAgICAgLy8gZWRnZXMuc2V0KHByZXZQb2ludCwgbmJQb2ludHMpO1xuICAgICAgICAgICAgICAgICAgICBlZGdlc1twcmV2UG9pbnRdID0gbmJQb2ludHM7XG4gICAgICAgICAgICAgICAgfVxuICAgICAgICAgICAgICAgIHByZXZQb2ludCA9IG5iUG9pbnRzO1xuICAgICAgICAgICAgICAgIG5iUG9pbnRzKys7XG4gICAgICAgICAgICAgICAgYnJlYWs7XG4gICAgICAgICAgICBjYXNlIFwielwiOlxuICAgICAgICAgICAgICAgIC8vIGVkZ2VzLnNldChwcmV2UG9pbnQsIGJlZ2luaW5nUGF0aCk7XG4gICAgICAgICAgICAgICAgZWRnZXNbcHJldlBvaW50XSA9IG5iUG9pbnRzO1xuICAgICAgICAgICAgICAgIGJlZ2luaW5nUGF0aCA9IHVuZGVmaW5lZDtcbiAgICAgICAgICAgICAgICBwcmV2UG9pbnQgPSB1bmRlZmluZWQ7XG4gICAgICAgICAgICAgICAgYnJlYWs7ICAgIFxuICAgICAgICB9XG4gICAgfVxuICAgIHJldHVybiB7cG9pbnRzIDogcG9pbnRzLCBlZGdlcyA6IGVkZ2VzfTtcbn1cblxuLy8gaW5pdGlhbGl6ZSBwb2ludHNcbnZhciBwb2ludHMgPSBbXTtcbnZhciBmb3JjZWRFZGdlcztcbnZhciBjaXR5U2V0O1xuXG5pZiAodGV4dE1lc2gpe1xuXG4gICAgdmFyIG15VGV4dCA9IHN2Z1RvUG9pbnRzKHN2Z1N0cmluZyk7XG4gICAgcG9pbnRzID0gbXlUZXh0LnBvaW50cztcbiAgICBmb3JjZWRFZGdlcyA9IG15VGV4dC5lZGdlcztcbiAgICAvLyBjaXR5U2V0ID0gbmV3IFNldChyYW5nZSgwLCBwb2ludHMubGVuZ3RoKSk7XG4gICAgY2l0eVNldCA9IHJhbmdlKDAsIHBvaW50cy5sZW5ndGgpO1xuXG4gICAgdmFyIHNjYWxlWCA9IDAuNTtcbiAgICB2YXIgc2NhbGVZID0gMC4yO1xuICAgIHZhciBkZWx0YVggPSAwLjI1O1xuICAgIHZhciBkZWx0YVkgPSAwLjM7XG5cbiAgICAvLyBzY2FsZSBwb2ludHMgdG8gWzAsMV0gKyBkZWx0YVxuICAgIHZhciBtYXhYID0gTWF0aC5tYXguYXBwbHkoTWF0aCwgcG9pbnRzLm1hcChmdW5jdGlvbihwKXtyZXR1cm4gcC54fSkpO1xuICAgIHZhciBtaW5YID0gTWF0aC5taW4uYXBwbHkoTWF0aCwgcG9pbnRzLm1hcChmdW5jdGlvbihwKXtyZXR1cm4gcC54fSkpO1xuICAgIHZhciBtYXhZID0gTWF0aC5tYXguYXBwbHkoTWF0aCwgcG9pbnRzLm1hcChmdW5jdGlvbihwKXtyZXR1cm4gcC55fSkpO1xuICAgIHZhciBtaW5ZID0gTWF0aC5taW4uYXBwbHkoTWF0aCwgcG9pbnRzLm1hcChmdW5jdGlvbihwKXtyZXR1cm4gcC55fSkpO1xuICAgIHBvaW50cyA9IHBvaW50cy5tYXAoZnVuY3Rpb24ocCl7XG4gICAgICAgIHZhciB4ID0gc2NhbGVYICogKHAueC1taW5YKS8obWF4WC1taW5YKSArIGRlbHRhWDtcbiAgICAgICAgdmFyIHkgPSBzY2FsZVkgKiAocC55LW1pblkpLyhtYXhZLW1pblkpICsgZGVsdGFZO1xuICAgICAgICB2YXIgbmV3UG9pbnQgPSBuZXcgUG9pbnQoeCwgeSk7XG4gICAgICAgIG5ld1BvaW50LmlkID0gcC5pZDtcblxuICAgICAgICByZXR1cm4gbmV3UG9pbnQ7XG4gICAgfSk7XG5cbiAgICAvLyBvbmx5IGFkZCByYW5kb20gcG9pbnRzXG4gICAgdmFyIG5iUG9pbnRzID0gcG9pbnRzLmxlbmd0aDtcbiAgICBmb3IodmFyIGk9MDsgaTxuYlJhbmRvbVBvaW50czsgKytpKSB7XG5cbiAgICAgICAgLy92YXIgeCA9IHJhbmRvbSgpICogdyArIHhJbml0O1xuICAgICAgICAvL3ZhciB5ID0gcmFuZG9tKCkgKiBoICsgeUluaXQ7XG5cbiAgICAgICAgdmFyIHggPSByYW5kb20oKTtcbiAgICAgICAgdmFyIHkgPSByYW5kb20oKTtcblxuICAgICAgICB2YXIgbmV3UG9pbnQgPSBuZXcgUG9pbnQoeCwgeSk7XG4gICAgICAgIG5ld1BvaW50LmlkID0gbmJQb2ludHM7XG5cbiAgICAgICAgcG9pbnRzLnB1c2gobmV3UG9pbnQpO1xuXG4gICAgICAgIG5iUG9pbnRzKys7XG4gICAgfVxuXG59IGVsc2Uge1xuICAgIC8vYWRkIHJhbmRvbSBwb2ludHNcblxuICAgIHZhciBuYlBvaW50cyA9IDA7XG4gICAgZm9yKHZhciBpPTA7IGk8bmJSYW5kb21Qb2ludHM7ICsraSkge1xuXG4gICAgICAgIHZhciB4ID0gcmFuZG9tKCk7XG4gICAgICAgIHZhciB5ID0gcmFuZG9tKCk7XG5cbiAgICAgICAgdmFyIG5ld1BvaW50ID0gbmV3IFBvaW50KHgsIHkpO1xuICAgICAgICBuZXdQb2ludC5pZCA9IG5iUG9pbnRzO1xuXG4gICAgICAgIHBvaW50cy5wdXNoKG5ld1BvaW50KTtcbiAgICAgICAgXG4gICAgICAgIG5iUG9pbnRzKys7XG4gICAgfVxuXG4gICAgLy8gY2l0eVNldCA9IG5ldyBTZXQocmFuZ2UoMCwgbmJDaXR5KSk7XG4gICAgY2l0eVNldCA9IHJhbmdlKDAsIG5iQ2l0eSk7XG4gICAgY29uc29sZS5sb2coY2l0eVNldCk7XG59XG5cblxuLy8gaW5pdGlhbGl6ZSBzdGFydCBwb2ludHNcbnZhciBwb3NzaWJsZVN0YXJ0UG9pbnRzSWQgPSBbXTtcblxuZm9yICh2YXIgaSA9IDA7IGkgPCBuYlN0YXJ0UG9pbnRzOyBpKyspe1xuICAgIHBvc3NpYmxlU3RhcnRQb2ludHNJZC5wdXNoKE1hdGguZmxvb3IobmJSYW5kb21Qb2ludHMgKiByYW5kb20oKSkpO1xufVxuXG5cblxubW9kdWxlLmV4cG9ydHMgPSB7XG4gICAgdGV4dE1lc2g6IHRleHRNZXNoLFxuICAgIHBvaW50czogcG9pbnRzLFxuICAgIGNpdHlTZXQ6IGNpdHlTZXQsXG4gICAgcG9zc2libGVTdGFydFBvaW50c0lkOiBwb3NzaWJsZVN0YXJ0UG9pbnRzSWQsXG4gICAgbmJSYW5kb21Qb2ludHM6IG5iUmFuZG9tUG9pbnRzLFxuICAgIGZvcmNlZEVkZ2VzOiBmb3JjZWRFZGdlc1xufSIsIid1c2Ugc3RyaWN0J1xuXG5tb2R1bGUuZXhwb3J0cyA9IGZ1bmN0aW9uIChjb250YWluZXIpe1xuXG5cdHZhciBtb3VzZSA9IHtcblx0ICAgIHg6IDAsXG5cdCAgICB5OiAwXG5cdH07XG5cblx0Y29udGFpbmVyLmFkZEV2ZW50TGlzdGVuZXIoICdtb3VzZW1vdmUnLCBmdW5jdGlvbihlKXtcblx0ICAgIHZhciByZWN0ID0gY29udGFpbmVyLmdldEJvdW5kaW5nQ2xpZW50UmVjdCgpO1xuXHQgICAgbW91c2UueCA9IGUuY2xpZW50WCAvIHJlY3Qud2lkdGg7XG5cdCAgICBtb3VzZS55ID0gZS5jbGllbnRZIC8gcmVjdC5oZWlnaHQ7XG5cdH0pO1xuXG5cdHJldHVybiBtb3VzZTtcblxufTtcbiIsIid1c2Ugc3RyaWN0J1xuXG5mdW5jdGlvbiBQb2ludCh4LCB5KSB7XG4gICAgdGhpcy5pZCA9IHVuZGVmaW5lZDsgICAgICAgICAgICAgICAgXG4gICAgdGhpcy54ID0geDtcbiAgICB0aGlzLnkgPSB5O1xuICAgIHRoaXMubmV4dHMgPSBbXTtcbn1cblxubW9kdWxlLmV4cG9ydHMgPSBQb2ludDsiLCIndXNlIHN0cmljdCdcblxudmFyIHJhbmRvbSA9IE1hdGgucmFuZG9tO1xuXG5tb2R1bGUuZXhwb3J0cyA9IGZ1bmN0aW9uKGNvbnRhaW5lcil7XG4gICAgXG4gICAgaWYoIWNvbnRhaW5lcilcbiAgICAgICAgdGhyb3cgbmV3IFR5cGVFcnJvcignTWlzc2luZyBjb250YWluZXInKTtcblxuICAgIHZhciBwb2ludHMgPSByZXF1aXJlKCcuL2luaXRpYWxpemVQb2ludHMuanMnKS5wb2ludHM7XG4gICAgdmFyIGNpdHlTZXQgPSByZXF1aXJlKCcuL2luaXRpYWxpemVQb2ludHMuanMnKS5jaXR5U2V0O1xuXG4gICAgdmFyIGVkZ2VzID0gcmVxdWlyZSgnLi9jcmVhdGVFZGdlcy5qcycpO1xuXG4gICAgdmFyIHBvcHVsYXRpb24gPSByZXF1aXJlKCcuL2luaXRpYWxpemVBbnRzJykoY29udGFpbmVyKTtcbiAgICB2YXIgbmJBbnRzID0gcG9wdWxhdGlvbi5sZW5ndGg7XG4gICAgICAgIFxuICAgIHZhciBjYW52YXMgPSBkb2N1bWVudC5jcmVhdGVFbGVtZW50KFwiY2FudmFzXCIpO1xuICAgIHZhciByZWN0ID0gY29udGFpbmVyLmdldEJvdW5kaW5nQ2xpZW50UmVjdCgpO1xuICAgIGNhbnZhcy53aWR0aCA9IHJlY3Qud2lkdGg7XG4gICAgY2FudmFzLmhlaWdodCA9IHJlY3QuaGVpZ2h0O1xuICAgIGNvbnRhaW5lci5hcHBlbmRDaGlsZChjYW52YXMpO1xuICAgIFxuICAgIHZhciBjb250ZXh0ID0gY2FudmFzLmdldENvbnRleHQoXCIyZFwiKTtcbiAgICBcblxuICAgIGZ1bmN0aW9uIHRpY2soKSB7XG4gICAgICAgIHZhciB3ID0gY2FudmFzLndpZHRoO1xuICAgICAgICB2YXIgaCA9IGNhbnZhcy5oZWlnaHQ7XG4gICAgICAgIHZhciBtb3VzZSA9IFtsYXN0TW91c2VNb3ZlRXZlbnQuY2xpZW50WC93LCBsYXN0TW91c2VNb3ZlRXZlbnQuY2xpZW50WS9oXTtcbiAgICAgICAgY29udGV4dC5zZXRUcmFuc2Zvcm0odywgMCwgMCwgaCwgMCwgMCk7XG4gICAgICAgIGNvbnRleHQuZmlsbFN0eWxlID0gXCIjZmZmXCI7XG4gICAgICAgIGNvbnRleHQuZmlsbFJlY3QoMCwwLHcsaCk7XG5cbiAgICAgICAgLy8gZWRnZXNcbiAgICAgICAgLy8gY29udGV4dC5zdHJva2VTdHlsZSA9IFwiIzAwMFwiO1xuICAgICAgICAvLyBmb3IodmFyIGk9MDsgaSA8IGVkZ2VzLmxlbmd0aDsgKytpKSB7XG4gICAgICAgIC8vICAgICBjb250ZXh0LmxpbmVXaWR0aCA9IDAuMDAwMTtcbiAgICAgICAgLy8gICAgIHZhciBlZGdlID0gZWRnZXNbaV07XG4gICAgICAgIC8vICAgICAvLyBpZiAoZWRnZS5waGVyb21vbiAhPSAwKXtcbiAgICAgICAgLy8gICAgIC8vICAgICBjb250ZXh0LmxpbmVXaWR0aCA9IE1hdGgubWluKDAuMDAwMDEgKiBlZGdlLnBoZXJvbW9uLCAwLjAxKTtcbiAgICAgICAgLy8gICAgIC8vIH0gZWxzZSB7XG4gICAgICAgIC8vICAgICAvLyAgICAgY29udGV4dC5saW5lV2lkdGggPSAwLjAwMDAxO1xuICAgICAgICAvLyAgICAgLy8gfVxuICAgICAgICAvLyAgICAgY29udGV4dC5iZWdpblBhdGgoKTtcbiAgICAgICAgLy8gICAgIGNvbnRleHQubW92ZVRvKHBvaW50c1tlZGdlLnB0MS5pZF0ueCwgcG9pbnRzW2VkZ2UucHQxLmlkXS55KTtcbiAgICAgICAgLy8gICAgIGNvbnRleHQubGluZVRvKHBvaW50c1tlZGdlLnB0Mi5pZF0ueCwgcG9pbnRzW2VkZ2UucHQyLmlkXS55KTtcbiAgICAgICAgLy8gICAgIGNvbnRleHQuc3Ryb2tlKCk7XG4gICAgICAgIC8vIH1cblxuICAgICAgICAvLyAvLyB2ZXJ0aWNlc1xuICAgICAgICAvLyBmb3IodmFyIGk9MDsgaTxwb2ludHMubGVuZ3RoOyArK2kpIHtcbiAgICAgICAgLy8gICAgIGNvbnRleHQuYmVnaW5QYXRoKClcbiAgICAgICAgLy8gICAgIHZhciBwb2ludCA9IHBvaW50c1tpXTtcbiAgICAgICAgLy8gICAgIGlmIChjaXR5U2V0Lmhhcyhwb2ludC5pZCkpIHtcbiAgICAgICAgLy8gICAgICAgICBjb250ZXh0LmZpbGxTdHlsZSA9IFwiIzAxMDFERlwiO1xuICAgICAgICAvLyAgICAgICAgIGNvbnRleHQuYXJjKHBvaW50LngsIHBvaW50LnksIDAuMDA2LCAwLCAyKk1hdGguUEkpO1xuICAgICAgICAvLyAgICAgfVxuICAgICAgICAvLyAgICAgZWxzZSB7XG4gICAgICAgIC8vICAgICAgICAgY29udGV4dC5maWxsU3R5bGUgPSBcIiMwMDBcIjtcbiAgICAgICAgLy8gICAgICAgICBjb250ZXh0LmFyYyhwb2ludHNbaV0ueCwgcG9pbnRzW2ldLnksIDAuMDAzLCAwLCAyKk1hdGguUEkpO1xuICAgICAgICAvLyAgICAgfVxuICAgICAgICAvLyAgICAgY29udGV4dC5jbG9zZVBhdGgoKTtcbiAgICAgICAgLy8gICAgIGNvbnRleHQuZmlsbCgpO1xuICAgICAgICAvLyB9XG5cbiAgICAgICAgLy8gbW92ZSBhbnRzXG4gICAgICAgIGZvciAoaSA9IDA7IGkgPCBuYkFudHM7IGkrKykge1xuICAgICAgICAgICAgcG9wdWxhdGlvbltpXS50cmFuc2l0KCk7XG4gICAgICAgIH1cblxuICAgICAgICAvLyBwaGVyb21vbiBldmFwb3JhdGlvblxuICAgICAgICBmb3IgKGkgPSAwOyBpIDwgZWRnZXMubGVuZ3RoOyBpKyspIHtcbiAgICAgICAgICAgIGlmKGVkZ2VzW2ldLnBoZXJvbW9uID4gMCl7XG4gICAgICAgICAgICAgICAgZWRnZXNbaV0ucGhlcm9tb24gLT0gMC4wMDAxO1xuICAgICAgICAgICAgfVxuICAgICAgICB9XG5cbiAgICAgICAgLy8gYW50c1xuICAgICAgICBmb3IodmFyIGk9MDsgaTxwb3B1bGF0aW9uLmxlbmd0aDsgKytpKSB7XG4gICAgICAgICAgICBjb250ZXh0LmJlZ2luUGF0aCgpXG4gICAgICAgICAgICB2YXIgeCA9IHBvcHVsYXRpb25baV0ueCArIDAuMDA1KnJhbmRvbSgpO1xuICAgICAgICAgICAgdmFyIHkgPSBwb3B1bGF0aW9uW2ldLnkgKyAwLjAwNSpyYW5kb20oKTtcbiAgICAgICAgICAgIC8vIHZhciB4ID0gcG9wdWxhdGlvbltpXS54O1xuICAgICAgICAgICAgLy8gdmFyIHkgPSBwb3B1bGF0aW9uW2ldLnk7XG4gICAgICAgICAgICBjb250ZXh0LmZpbGxTdHlsZSA9IFwiYmxhY2tcIlxuICAgICAgICAgICAgY29udGV4dC5maWxsUmVjdCh4LCB5LCAwLjAwMTIsIDAuMDAxMik7XG4gICAgICAgICAgICBjb250ZXh0LmNsb3NlUGF0aCgpO1xuICAgICAgICAgICAgY29udGV4dC5maWxsKCk7XG4gICAgICAgIH1cblxuICAgIH07XG4gICAgXG4gICAgdmFyIGxhc3RNb3VzZU1vdmVFdmVudCA9IHtcbiAgICAgICAgY2xpZW50WDogMCxcbiAgICAgICAgY2xpZW50WTogMFxuICAgIH07XG4gICAgXG4gICAgY29udGFpbmVyLmFkZEV2ZW50TGlzdGVuZXIoJ21vdXNlbW92ZScsIGZ1bmN0aW9uKGUpe1xuICAgICAgICBsYXN0TW91c2VNb3ZlRXZlbnQgPSBlO1xuICAgIH0pO1xuICAgIFxuICAgIHZhciBwYXVzZWQgPSBmYWxzZTtcbiAgICBcbiAgICBmdW5jdGlvbiB0b2dnbGVQbGF5UGF1c2UoKXtcbiAgICAgICAgcGF1c2VkID0gIXBhdXNlZDtcbiAgICAgICAgaWYoIXBhdXNlZClcbiAgICAgICAgICAgIGFuaW1hdGUoKTtcbiAgICB9XG4gICAgXG4gICAgY29udGFpbmVyLmFkZEV2ZW50TGlzdGVuZXIoJ2NsaWNrJywgdG9nZ2xlUGxheVBhdXNlKTtcbiAgICBcbiAgICBmdW5jdGlvbiBhbmltYXRlKCl7XG4gICAgICAgIHRpY2soKTtcbiAgICAgICAgXG4gICAgICAgIGlmKCFwYXVzZWQpXG4gICAgICAgICAgICByZXF1ZXN0QW5pbWF0aW9uRnJhbWUoYW5pbWF0ZSk7XG4gICAgfTtcbiAgICBhbmltYXRlKCk7XG4gICAgXG4gICAgcmV0dXJuIHtcbiAgICAgICAgdG9nZ2xlUGxheVBhdXNlOiB0b2dnbGVQbGF5UGF1c2UsXG4gICAgICAgIC8vIHNob3VsZCBiZSBhIGdldHRlci9zZXR0ZXIsIGJ1dCBJRThcbiAgICAgICAgZ2V0QW50Q291bnQ6IGZ1bmN0aW9uKCl7XG4gICAgICAgICAgICB0aHJvdyAnVE9ETyc7XG4gICAgICAgIH0sXG4gICAgICAgIHNldEFudENvdW50OiBmdW5jdGlvbigpe1xuICAgICAgICAgICAgdGhyb3cgJ1RPRE8nO1xuICAgICAgICB9XG4gICAgfVxufVxuIiwiJ3VzZSBzdHJpY3QnO1xuXG5yZXF1aXJlKCcuLi9pbmRleC5qcycpKGRvY3VtZW50LmJvZHkpOyIsIid1c2Ugc3RyaWN0J1xuXG52YXIgc3FydCA9IE1hdGguc3FydDtcbnZhciBwb3cgPSBNYXRoLnBvdztcblxuZnVuY3Rpb24gc2lnbih4KSB7XG5cdHJldHVybiB4ID8geCA8IDAgPyAtMSA6IDEgOiAwO1xufVxuXG5mdW5jdGlvbiByYW5nZShzdGFydCwgY291bnQpIHtcbiAgICByZXR1cm4gQXJyYXkuYXBwbHkoMCwgQXJyYXkoY291bnQpKS5tYXAoZnVuY3Rpb24gKGVsZW1lbnQsIGluZGV4KSB7XG4gICAgXHRyZXR1cm4gaW5kZXggKyBzdGFydFxuICAgIH0pO1xufVxuXG5mdW5jdGlvbiBkaXN0YW5jZShhLCBiKXtcblx0cmV0dXJuIHNxcnQocG93KGEueCAtIGIueCwgMikgKyBwb3coYS55IC0gYi55LCAyKSk7XG59XG5cbmZ1bmN0aW9uIG5vcm0odil7XG5cdHJldHVybiBzcXJ0KHBvdyh2LngsIDIpICsgcG93KHYueSwgMikpO1xufVxuXG5tb2R1bGUuZXhwb3J0cyA9IHtcblx0c2lnbjogc2lnbixcblx0cmFuZ2U6IHJhbmdlLFxuXHRkaXN0YW5jZTogZGlzdGFuY2UsXG5cdG5vcm06IG5vcm1cbn0iLCIndXNlIHN0cmljdCdcblxuZnVuY3Rpb24gVmVjdG9yKHgsIHkpIHtcbiAgICB0aGlzLnggPSB4OyAgICAgICAgICAgICAgICBcbiAgICB0aGlzLnkgPSB5O1xufVxuXG5WZWN0b3IucHJvdG90eXBlLm5vcm0gPSBmdW5jdGlvbigpe1xuXHRyZXR1cm4gTWF0aC5zcXJ0KHRoaXMueCAqIHRoaXMueCArIHRoaXMueSAqIHRoaXMueSk7XG59XG5cblZlY3Rvci5wcm90b3R5cGUubm9ybWFsaXplID0gZnVuY3Rpb24oKXtcblx0dmFyIG5vcm0gPSB0aGlzLm5vcm0oKTtcblx0dGhpcy54ID0gdGhpcy54IC8gbm9ybTtcblx0dGhpcy55ID0gdGhpcy55IC8gbm9ybTtcbn1cblxuXG5cbm1vZHVsZS5leHBvcnRzID0gVmVjdG9yOyJdfQ==