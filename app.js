(function e(t,n,r){function s(o,u){if(!n[o]){if(!t[o]){var a=typeof require=="function"&&require;if(!u&&a)return a(o,!0);if(i)return i(o,!0);var f=new Error("Cannot find module '"+o+"'");throw f.code="MODULE_NOT_FOUND",f}var l=n[o]={exports:{}};t[o][0].call(l.exports,function(e){var n=t[o][1][e];return s(n?n:e)},l,l.exports,e,t,n,r)}return n[o].exports}var i=typeof require=="function"&&require;for(var o=0;o<r.length;o++)s(r[o]);return s})({"/Users/Romain/Documents/Programmation/Projets_Web/AntColony/node_modules/delaunay-triangulate/node_modules/incremental-convex-hull/ich.js":[function(require,module,exports){
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
},{"robust-orientation":"/Users/Romain/Documents/Programmation/Projets_Web/AntColony/node_modules/delaunay-triangulate/node_modules/incremental-convex-hull/node_modules/robust-orientation/orientation.js","simplicial-complex":"/Users/Romain/Documents/Programmation/Projets_Web/AntColony/node_modules/delaunay-triangulate/node_modules/incremental-convex-hull/node_modules/simplicial-complex/topology.js"}],"/Users/Romain/Documents/Programmation/Projets_Web/AntColony/node_modules/delaunay-triangulate/node_modules/incremental-convex-hull/node_modules/robust-orientation/node_modules/robust-scale/node_modules/two-sum/two-sum.js":[function(require,module,exports){
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
},{}],"/Users/Romain/Documents/Programmation/Projets_Web/AntColony/node_modules/delaunay-triangulate/node_modules/incremental-convex-hull/node_modules/robust-orientation/node_modules/robust-scale/robust-scale.js":[function(require,module,exports){
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
},{"two-product":"/Users/Romain/Documents/Programmation/Projets_Web/AntColony/node_modules/delaunay-triangulate/node_modules/incremental-convex-hull/node_modules/robust-orientation/node_modules/two-product/two-product.js","two-sum":"/Users/Romain/Documents/Programmation/Projets_Web/AntColony/node_modules/delaunay-triangulate/node_modules/incremental-convex-hull/node_modules/robust-orientation/node_modules/robust-scale/node_modules/two-sum/two-sum.js"}],"/Users/Romain/Documents/Programmation/Projets_Web/AntColony/node_modules/delaunay-triangulate/node_modules/incremental-convex-hull/node_modules/robust-orientation/node_modules/robust-subtract/robust-diff.js":[function(require,module,exports){
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
},{}],"/Users/Romain/Documents/Programmation/Projets_Web/AntColony/node_modules/delaunay-triangulate/node_modules/incremental-convex-hull/node_modules/robust-orientation/node_modules/robust-sum/robust-sum.js":[function(require,module,exports){
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
},{}],"/Users/Romain/Documents/Programmation/Projets_Web/AntColony/node_modules/delaunay-triangulate/node_modules/incremental-convex-hull/node_modules/robust-orientation/node_modules/two-product/two-product.js":[function(require,module,exports){
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
},{}],"/Users/Romain/Documents/Programmation/Projets_Web/AntColony/node_modules/delaunay-triangulate/node_modules/incremental-convex-hull/node_modules/robust-orientation/orientation.js":[function(require,module,exports){
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
},{"robust-scale":"/Users/Romain/Documents/Programmation/Projets_Web/AntColony/node_modules/delaunay-triangulate/node_modules/incremental-convex-hull/node_modules/robust-orientation/node_modules/robust-scale/robust-scale.js","robust-subtract":"/Users/Romain/Documents/Programmation/Projets_Web/AntColony/node_modules/delaunay-triangulate/node_modules/incremental-convex-hull/node_modules/robust-orientation/node_modules/robust-subtract/robust-diff.js","robust-sum":"/Users/Romain/Documents/Programmation/Projets_Web/AntColony/node_modules/delaunay-triangulate/node_modules/incremental-convex-hull/node_modules/robust-orientation/node_modules/robust-sum/robust-sum.js","two-product":"/Users/Romain/Documents/Programmation/Projets_Web/AntColony/node_modules/delaunay-triangulate/node_modules/incremental-convex-hull/node_modules/robust-orientation/node_modules/two-product/two-product.js"}],"/Users/Romain/Documents/Programmation/Projets_Web/AntColony/node_modules/delaunay-triangulate/node_modules/incremental-convex-hull/node_modules/simplicial-complex/node_modules/bit-twiddle/twiddle.js":[function(require,module,exports){
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


},{}],"/Users/Romain/Documents/Programmation/Projets_Web/AntColony/node_modules/delaunay-triangulate/node_modules/incremental-convex-hull/node_modules/simplicial-complex/node_modules/union-find/index.js":[function(require,module,exports){
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
},{}],"/Users/Romain/Documents/Programmation/Projets_Web/AntColony/node_modules/delaunay-triangulate/node_modules/incremental-convex-hull/node_modules/simplicial-complex/topology.js":[function(require,module,exports){
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

},{"bit-twiddle":"/Users/Romain/Documents/Programmation/Projets_Web/AntColony/node_modules/delaunay-triangulate/node_modules/incremental-convex-hull/node_modules/simplicial-complex/node_modules/bit-twiddle/twiddle.js","union-find":"/Users/Romain/Documents/Programmation/Projets_Web/AntColony/node_modules/delaunay-triangulate/node_modules/incremental-convex-hull/node_modules/simplicial-complex/node_modules/union-find/index.js"}],"/Users/Romain/Documents/Programmation/Projets_Web/AntColony/node_modules/delaunay-triangulate/node_modules/uniq/uniq.js":[function(require,module,exports){
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

},{}],"/Users/Romain/Documents/Programmation/Projets_Web/AntColony/node_modules/delaunay-triangulate/triangulate.js":[function(require,module,exports){
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
},{"incremental-convex-hull":"/Users/Romain/Documents/Programmation/Projets_Web/AntColony/node_modules/delaunay-triangulate/node_modules/incremental-convex-hull/ich.js","uniq":"/Users/Romain/Documents/Programmation/Projets_Web/AntColony/node_modules/delaunay-triangulate/node_modules/uniq/uniq.js"}],"/Users/Romain/Documents/Programmation/Projets_Web/AntColony/node_modules/game-shell/lib/hrtime-polyfill.js":[function(require,module,exports){
if(typeof window.performance === "object") {
  if(window.performance.now) {
    module.exports = function() { return window.performance.now() }
  } else if(window.performance.webkitNow) {
    module.exports = function() { return window.performance.webkitNow() }
  }
} else if(Date.now) {
  module.exports = Date.now
} else {
  module.exports = function() { return (new Date()).getTime() }
}

},{}],"/Users/Romain/Documents/Programmation/Projets_Web/AntColony/node_modules/game-shell/lib/mousewheel-polyfill.js":[function(require,module,exports){
//Adapted from here: https://developer.mozilla.org/en-US/docs/Web/Reference/Events/wheel?redirectlocale=en-US&redirectslug=DOM%2FMozilla_event_reference%2Fwheel

var prefix = "", _addEventListener, onwheel, support;

// detect event model
if ( window.addEventListener ) {
  _addEventListener = "addEventListener";
} else {
  _addEventListener = "attachEvent";
  prefix = "on";
}

// detect available wheel event
support = "onwheel" in document.createElement("div") ? "wheel" : // Modern browsers support "wheel"
          document.onmousewheel !== undefined ? "mousewheel" : // Webkit and IE support at least "mousewheel"
          "DOMMouseScroll"; // let's assume that remaining browsers are older Firefox

function _addWheelListener( elem, eventName, callback, useCapture ) {
  elem[ _addEventListener ]( prefix + eventName, support == "wheel" ? callback : function( originalEvent ) {
    !originalEvent && ( originalEvent = window.event );

    // create a normalized event object
    var event = {
      // keep a ref to the original event object
      originalEvent: originalEvent,
      target: originalEvent.target || originalEvent.srcElement,
      type: "wheel",
      deltaMode: originalEvent.type == "MozMousePixelScroll" ? 0 : 1,
      deltaX: 0,
      delatZ: 0,
      preventDefault: function() {
        originalEvent.preventDefault ?
          originalEvent.preventDefault() :
          originalEvent.returnValue = false;
      }
    };
    
    // calculate deltaY (and deltaX) according to the event
    if ( support == "mousewheel" ) {
      event.deltaY = - 1/40 * originalEvent.wheelDelta;
      // Webkit also support wheelDeltaX
      originalEvent.wheelDeltaX && ( event.deltaX = - 1/40 * originalEvent.wheelDeltaX );
    } else {
      event.deltaY = originalEvent.detail;
    }

    // it's time to fire the callback
    return callback( event );
  }, useCapture || false );
}

module.exports = function( elem, callback, useCapture ) {
  _addWheelListener( elem, support, callback, useCapture );

  // handle MozMousePixelScroll in older Firefox
  if( support == "DOMMouseScroll" ) {
    _addWheelListener( elem, "MozMousePixelScroll", callback, useCapture );
  }
};
},{}],"/Users/Romain/Documents/Programmation/Projets_Web/AntColony/node_modules/game-shell/lib/raf-polyfill.js":[function(require,module,exports){
// http://paulirish.com/2011/requestanimationframe-for-smart-animating/
// http://my.opera.com/emoller/blog/2011/12/20/requestanimationframe-for-smart-er-animating
 
// requestAnimationFrame polyfill by Erik Möller. fixes from Paul Irish and Tino Zijdel
 
// MIT license
var lastTime = 0;
var vendors = ['ms', 'moz', 'webkit', 'o'];
for(var x = 0; x < vendors.length && !window.requestAnimationFrame; ++x) {
    window.requestAnimationFrame = window[vendors[x]+'RequestAnimationFrame'];
    window.cancelAnimationFrame = window[vendors[x]+'CancelAnimationFrame'] 
                               || window[vendors[x]+'CancelRequestAnimationFrame'];
}

if (!window.requestAnimationFrame)
    window.requestAnimationFrame = function(callback, element) {
        var currTime = new Date().getTime();
        var timeToCall = Math.max(0, 16 - (currTime - lastTime));
        var id = window.setTimeout(function() { callback(currTime + timeToCall); }, 
          timeToCall);
        lastTime = currTime + timeToCall;
        return id;
    };

if (!window.cancelAnimationFrame)
    window.cancelAnimationFrame = function(id) {
        clearTimeout(id);
    };

},{}],"/Users/Romain/Documents/Programmation/Projets_Web/AntColony/node_modules/game-shell/node_modules/binary-search-bounds/search-bounds.js":[function(require,module,exports){
"use strict"

function compileSearch(funcName, predicate, reversed, extraArgs, useNdarray, earlyOut) {
  var code = [
    "function ", funcName, "(a,l,h,", extraArgs.join(","),  "){",
earlyOut ? "" : "var i=", (reversed ? "l-1" : "h+1"),
";while(l<=h){\
var m=(l+h)>>>1,x=a", useNdarray ? ".get(m)" : "[m]"]
  if(earlyOut) {
    if(predicate.indexOf("c") < 0) {
      code.push(";if(x===y){return m}else if(x<=y){")
    } else {
      code.push(";var p=c(x,y);if(p===0){return m}else if(p<=0){")
    }
  } else {
    code.push(";if(", predicate, "){i=m;")
  }
  if(reversed) {
    code.push("l=m+1}else{h=m-1}")
  } else {
    code.push("h=m-1}else{l=m+1}")
  }
  code.push("}")
  if(earlyOut) {
    code.push("return -1};")
  } else {
    code.push("return i};")
  }
  return code.join("")
}

function compileBoundsSearch(predicate, reversed, suffix, earlyOut) {
  var result = new Function([
  compileSearch("A", "x" + predicate + "y", reversed, ["y"], false, earlyOut),
  compileSearch("B", "x" + predicate + "y", reversed, ["y"], true, earlyOut),
  compileSearch("P", "c(x,y)" + predicate + "0", reversed, ["y", "c"], false, earlyOut),
  compileSearch("Q", "c(x,y)" + predicate + "0", reversed, ["y", "c"], true, earlyOut),
"function dispatchBsearch", suffix, "(a,y,c,l,h){\
if(a.shape){\
if(typeof(c)==='function'){\
return Q(a,(l===undefined)?0:l|0,(h===undefined)?a.shape[0]-1:h|0,y,c)\
}else{\
return B(a,(c===undefined)?0:c|0,(l===undefined)?a.shape[0]-1:l|0,y)\
}}else{\
if(typeof(c)==='function'){\
return P(a,(l===undefined)?0:l|0,(h===undefined)?a.length-1:h|0,y,c)\
}else{\
return A(a,(c===undefined)?0:c|0,(l===undefined)?a.length-1:l|0,y)\
}}}\
return dispatchBsearch", suffix].join(""))
  return result()
}

module.exports = {
  ge: compileBoundsSearch(">=", false, "GE"),
  gt: compileBoundsSearch(">", false, "GT"),
  lt: compileBoundsSearch("<", true, "LT"),
  le: compileBoundsSearch("<=", true, "LE"),
  eq: compileBoundsSearch("-", true, "EQ", true)
}

},{}],"/Users/Romain/Documents/Programmation/Projets_Web/AntColony/node_modules/game-shell/node_modules/domready/ready.js":[function(require,module,exports){
/*!
  * domready (c) Dustin Diaz 2014 - License MIT
  */
!function (name, definition) {

  if (typeof module != 'undefined') module.exports = definition()
  else if (typeof define == 'function' && typeof define.amd == 'object') define(definition)
  else this[name] = definition()

}('domready', function () {

  var fns = [], listener
    , doc = document
    , domContentLoaded = 'DOMContentLoaded'
    , loaded = /^loaded|^c/.test(doc.readyState)

  if (!loaded)
  doc.addEventListener(domContentLoaded, listener = function () {
    doc.removeEventListener(domContentLoaded, listener)
    loaded = 1
    while (listener = fns.shift()) listener()
  })

  return function (fn) {
    loaded ? fn() : fns.push(fn)
  }

});

},{}],"/Users/Romain/Documents/Programmation/Projets_Web/AntColony/node_modules/game-shell/node_modules/invert-hash/invert.js":[function(require,module,exports){
"use strict"

function invert(hash) {
  var result = {}
  for(var i in hash) {
    if(hash.hasOwnProperty(i)) {
      result[hash[i]] = i
    }
  }
  return result
}

module.exports = invert
},{}],"/Users/Romain/Documents/Programmation/Projets_Web/AntColony/node_modules/game-shell/node_modules/iota-array/iota.js":[function(require,module,exports){
"use strict"

function iota(n) {
  var result = new Array(n)
  for(var i=0; i<n; ++i) {
    result[i] = i
  }
  return result
}

module.exports = iota
},{}],"/Users/Romain/Documents/Programmation/Projets_Web/AntColony/node_modules/game-shell/node_modules/uniq/uniq.js":[function(require,module,exports){
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

},{}],"/Users/Romain/Documents/Programmation/Projets_Web/AntColony/node_modules/game-shell/node_modules/vkey/index.js":[function(require,module,exports){
var ua = typeof window !== 'undefined' ? window.navigator.userAgent : ''
  , isOSX = /OS X/.test(ua)
  , isOpera = /Opera/.test(ua)
  , maybeFirefox = !/like Gecko/.test(ua) && !isOpera

var i, output = module.exports = {
  0:  isOSX ? '<menu>' : '<UNK>'
, 1:  '<mouse 1>'
, 2:  '<mouse 2>'
, 3:  '<break>'
, 4:  '<mouse 3>'
, 5:  '<mouse 4>'
, 6:  '<mouse 5>'
, 8:  '<backspace>'
, 9:  '<tab>'
, 12: '<clear>'
, 13: '<enter>'
, 16: '<shift>'
, 17: '<control>'
, 18: '<alt>'
, 19: '<pause>'
, 20: '<caps-lock>'
, 21: '<ime-hangul>'
, 23: '<ime-junja>'
, 24: '<ime-final>'
, 25: '<ime-kanji>'
, 27: '<escape>'
, 28: '<ime-convert>'
, 29: '<ime-nonconvert>'
, 30: '<ime-accept>'
, 31: '<ime-mode-change>'
, 27: '<escape>'
, 32: '<space>'
, 33: '<page-up>'
, 34: '<page-down>'
, 35: '<end>'
, 36: '<home>'
, 37: '<left>'
, 38: '<up>'
, 39: '<right>'
, 40: '<down>'
, 41: '<select>'
, 42: '<print>'
, 43: '<execute>'
, 44: '<snapshot>'
, 45: '<insert>'
, 46: '<delete>'
, 47: '<help>'
, 91: '<meta>'  // meta-left -- no one handles left and right properly, so we coerce into one.
, 92: '<meta>'  // meta-right
, 93: isOSX ? '<meta>' : '<menu>'      // chrome,opera,safari all report this for meta-right (osx mbp).
, 95: '<sleep>'
, 106: '<num-*>'
, 107: '<num-+>'
, 108: '<num-enter>'
, 109: '<num-->'
, 110: '<num-.>'
, 111: '<num-/>'
, 144: '<num-lock>'
, 145: '<scroll-lock>'
, 160: '<shift-left>'
, 161: '<shift-right>'
, 162: '<control-left>'
, 163: '<control-right>'
, 164: '<alt-left>'
, 165: '<alt-right>'
, 166: '<browser-back>'
, 167: '<browser-forward>'
, 168: '<browser-refresh>'
, 169: '<browser-stop>'
, 170: '<browser-search>'
, 171: '<browser-favorites>'
, 172: '<browser-home>'

  // ff/osx reports '<volume-mute>' for '-'
, 173: isOSX && maybeFirefox ? '-' : '<volume-mute>'
, 174: '<volume-down>'
, 175: '<volume-up>'
, 176: '<next-track>'
, 177: '<prev-track>'
, 178: '<stop>'
, 179: '<play-pause>'
, 180: '<launch-mail>'
, 181: '<launch-media-select>'
, 182: '<launch-app 1>'
, 183: '<launch-app 2>'
, 186: ';'
, 187: '='
, 188: ','
, 189: '-'
, 190: '.'
, 191: '/'
, 192: '`'
, 219: '['
, 220: '\\'
, 221: ']'
, 222: "'"
, 223: '<meta>'
, 224: '<meta>'       // firefox reports meta here.
, 226: '<alt-gr>'
, 229: '<ime-process>'
, 231: isOpera ? '`' : '<unicode>'
, 246: '<attention>'
, 247: '<crsel>'
, 248: '<exsel>'
, 249: '<erase-eof>'
, 250: '<play>'
, 251: '<zoom>'
, 252: '<no-name>'
, 253: '<pa-1>'
, 254: '<clear>'
}

for(i = 58; i < 65; ++i) {
  output[i] = String.fromCharCode(i)
}

// 0-9
for(i = 48; i < 58; ++i) {
  output[i] = (i - 48)+''
}

// A-Z
for(i = 65; i < 91; ++i) {
  output[i] = String.fromCharCode(i)
}

// num0-9
for(i = 96; i < 106; ++i) {
  output[i] = '<num-'+(i - 96)+'>'
}

// F1-F24
for(i = 112; i < 136; ++i) {
  output[i] = 'F'+(i-111)
}

},{}],"/Users/Romain/Documents/Programmation/Projets_Web/AntColony/node_modules/game-shell/shell.js":[function(require,module,exports){
"use strict"

var EventEmitter = require("events").EventEmitter
  , util         = require("util")
  , domready     = require("domready")
  , vkey         = require("vkey")
  , invert       = require("invert-hash")
  , uniq         = require("uniq")
  , bsearch      = require("binary-search-bounds")
  , iota         = require("iota-array")
  , min          = Math.min

//Browser compatibility hacks
require("./lib/raf-polyfill.js")
var addMouseWheel = require("./lib/mousewheel-polyfill.js")
var hrtime = require("./lib/hrtime-polyfill.js")

//Remove angle braces and other useless crap
var filtered_vkey = (function() {
  var result = new Array(256)
    , i, j, k
  for(i=0; i<256; ++i) {
    result[i] = "UNK"
  }
  for(i in vkey) {
    k = vkey[i]
    if(k.charAt(0) === '<' && k.charAt(k.length-1) === '>') {
      k = k.substring(1, k.length-1)
    }
    k = k.replace(/\s/g, "-")
    result[parseInt(i)] = k
  }
  return result
})()

//Compute minimal common set of keyboard functions
var keyNames = uniq(Object.keys(invert(filtered_vkey)))

//Translates a virtual keycode to a normalized keycode
function virtualKeyCode(key) {
  return bsearch.eq(keyNames, key)
}

//Maps a physical keycode to a normalized keycode
function physicalKeyCode(key) {
  return virtualKeyCode(filtered_vkey[key])
}

//Game shell
function GameShell() {
  EventEmitter.call(this)
  this._curKeyState  = new Array(keyNames.length)
  this._pressCount   = new Array(keyNames.length)
  this._releaseCount = new Array(keyNames.length)
  
  this._tickInterval = null
  this._rafHandle = null
  this._tickRate = 0
  this._lastTick = hrtime()
  this._frameTime = 0.0
  this._paused = true
  this._width = 0
  this._height = 0
  
  this._wantFullscreen = false
  this._wantPointerLock = false
  this._fullscreenActive = false
  this._pointerLockActive = false
  
  this._render = render.bind(undefined, this)
  
  for(var i=0; i<keyNames.length; ++i) {
    this._curKeyState[i] = false
    this._pressCount[i] = this._releaseCount[i] = 0
  }
  
  //Public members
  this.element = null
  this.bindings = {}
  this.frameSkip = 100.0
  this.tickCount = 0
  this.frameCount = 0
  this.startTime = hrtime()
  this.tickTime = this._tickRate
  this.frameTime = 10.0
  this.stickyFullscreen = false
  this.stickyPointLock = false
  
  //Scroll stuff
  this.scroll = [0,0,0]
    
  //Mouse state
  this.mouseX = 0
  this.mouseY = 0
  this.prevMouseX = 0
  this.prevMouseY = 0
}

util.inherits(GameShell, EventEmitter)

var proto = GameShell.prototype

//Bind keynames
proto.keyNames = keyNames

//Binds a virtual keyboard event to a physical key
proto.bind = function(virtual_key) {
  //Look up previous key bindings
  var arr
  if(virtual_key in this.bindings) {
    arr = this.bindings[virtual_key]
  } else {
    arr = []
  }
  //Add keys to list
  var physical_key
  for(var i=1, n=arguments.length; i<n; ++i) {
    physical_key = arguments[i]
    if(virtualKeyCode(physical_key) >= 0) {
      arr.push(physical_key)
    } else if(physical_key in this.bindings) {
      var keybinds = this.bindings[physical_key]
      for(var j=0; j<keybinds.length; ++j) {
        arr.push(keybinds[j])
      }
    }
  }
  //Remove any duplicate keys
  arr = uniq(arr)
  if(arr.length > 0) {
    this.bindings[virtual_key] = arr
  }
  this.emit('bind', virtual_key, arr)
}

//Unbinds a virtual keyboard event
proto.unbind = function(virtual_key) {
  if(virtual_key in this.bindings) {
    delete this.bindings[virtual_key]
  }
  this.emit('unbind', virtual_key)
}

//Checks if a key is set in a given state
function lookupKey(state, bindings, key) {
  if(key in bindings) {
    var arr = bindings[key]
    for(var i=0, n=arr.length; i<n; ++i) {
      if(state[virtualKeyCode(arr[i])]) {
        return true
      }
    }
    return false
  }
  var kc = virtualKeyCode(key)
  if(kc >= 0) {
    return state[kc]
  }
  return false
}

//Checks if a key is set in a given state
function lookupCount(state, bindings, key) {
  if(key in bindings) {
    var arr = bindings[key], r = 0
    for(var i=0, n=arr.length; i<n; ++i) {
      r += state[virtualKeyCode(arr[i])]
    }
    return r
  }
  var kc = virtualKeyCode(key)
  if(kc >= 0) {
    return state[kc]
  }
  return 0
}

//Checks if a key (either physical or virtual) is currently held down
proto.down = function(key) {
  return lookupKey(this._curKeyState, this.bindings, key)
}

//Checks if a key was ever down
proto.wasDown = function(key) {
  return this.down(key) || !!this.press(key)
}

//Opposite of down
proto.up = function(key) {
  return !this.down(key)
}

//Checks if a key was released during previous frame
proto.wasUp = function(key) {
  return this.up(key) || !!this.release(key)
}

//Returns the number of times a key was pressed since last tick
proto.press = function(key) {
  return lookupCount(this._pressCount, this.bindings, key)
}

//Returns the number of times a key was released since last tick
proto.release = function(key) {
  return lookupCount(this._releaseCount, this.bindings, key)
}

//Pause/unpause the game loop
Object.defineProperty(proto, "paused", {
  get: function() {
    return this._paused
  },
  set: function(state) {
    var ns = !!state
    if(ns !== this._paused) {
      if(!this._paused) {
        this._paused = true
        this._frameTime = min(1.0, (hrtime() - this._lastTick) / this._tickRate)
        clearInterval(this._tickInterval)
        //cancelAnimationFrame(this._rafHandle)
      } else {
        this._paused = false
        this._lastTick = hrtime() - Math.floor(this._frameTime * this._tickRate)
        this._tickInterval = setInterval(tick, this._tickRate, this)
        this._rafHandle = requestAnimationFrame(this._render)
      }
    }
  }
})

//Fullscreen state toggle

function tryFullscreen(shell) {
  //Request full screen
  var elem = shell.element
  
  if(shell._wantFullscreen && !shell._fullscreenActive) {
    var fs = elem.requestFullscreen ||
             elem.requestFullScreen ||
             elem.webkitRequestFullscreen ||
             elem.webkitRequestFullScreen ||
             elem.mozRequestFullscreen ||
             elem.mozRequestFullScreen ||
             function() {}
    fs.call(elem)
  }
  if(shell._wantPointerLock && !shell._pointerLockActive) {
    var pl =  elem.requestPointerLock ||
              elem.webkitRequestPointerLock ||
              elem.mozRequestPointerLock ||
              elem.msRequestPointerLock ||
              elem.oRequestPointerLock ||
              function() {}
    pl.call(elem)
  }
}

var cancelFullscreen = document.exitFullscreen ||
                       document.cancelFullscreen ||  //Why can no one agree on this?
                       document.cancelFullScreen ||
                       document.webkitCancelFullscreen ||
                       document.webkitCancelFullScreen ||
                       document.mozCancelFullscreen ||
                       document.mozCancelFullScreen ||
                       function(){}

Object.defineProperty(proto, "fullscreen", {
  get: function() {
    return this._fullscreenActive
  },
  set: function(state) {
    var ns = !!state
    if(!ns) {
      this._wantFullscreen = false
      cancelFullscreen.call(document)
    } else {
      this._wantFullscreen = true
      tryFullscreen(this)
    }
    return this._fullscreenActive
  }
})

function handleFullscreen(shell) {
  shell._fullscreenActive = document.fullscreen ||
                            document.mozFullScreen ||
                            document.webkitIsFullScreen ||
                            false
  if(!shell.stickyFullscreen && shell._fullscreenActive) {
    shell._wantFullscreen = false
  }
}

//Pointer lock state toggle
var exitPointerLock = document.exitPointerLock ||
                      document.webkitExitPointerLock ||
                      document.mozExitPointerLock ||
                      function() {}

Object.defineProperty(proto, "pointerLock", {
  get: function() {
    return this._pointerLockActive
  },
  set: function(state) {
    var ns = !!state
    if(!ns) {
      this._wantPointerLock = false
      exitPointerLock.call(document)
    } else {
      this._wantPointerLock = true
      tryFullscreen(this)
    }
    return this._pointerLockActive
  }
})

function handlePointerLockChange(shell, event) {
  shell._pointerLockActive = shell.element === (
      document.pointerLockElement ||
      document.mozPointerLockElement ||
      document.webkitPointerLockElement ||
      null)
  if(!shell.stickyPointerLock && shell._pointerLockActive) {
    shell._wantPointerLock = false
  }
}

//Width and height
Object.defineProperty(proto, "width", {
  get: function() {
    return this.element.clientWidth
  }
})
Object.defineProperty(proto, "height", {
  get: function() {
    return this.element.clientHeight
  }
})

//Set key state
function setKeyState(shell, key, state) {
  var ps = shell._curKeyState[key]
  if(ps !== state) {
    if(state) {
      shell._pressCount[key]++
    } else {
      shell._releaseCount[key]++
    }
    shell._curKeyState[key] = state
  }
}

//Ticks the game state one update
function tick(shell) {
  var skip = hrtime() + shell.frameSkip
    , pCount = shell._pressCount
    , rCount = shell._releaseCount
    , i, s, t
    , tr = shell._tickRate
    , n = keyNames.length
  while(!shell._paused &&
        hrtime() >= shell._lastTick + tr) {
    
    //Skip frames if we are over budget
    if(hrtime() > skip) {
      shell._lastTick = hrtime() + tr
      return
    }
    
    //Tick the game
    s = hrtime()
    shell.emit("tick")
    t = hrtime()
    shell.tickTime = t - s
    
    //Update counters and time
    ++shell.tickCount
    shell._lastTick += tr
    
    //Shift input state
    for(i=0; i<n; ++i) {
      pCount[i] = rCount[i] = 0
    }
    if(shell._pointerLockActive) {
      shell.prevMouseX = shell.mouseX = shell.width>>1
      shell.prevMouseY = shell.mouseY = shell.height>>1
    } else {
      shell.prevMouseX = shell.mouseX
      shell.prevMouseY = shell.mouseY
    }
    shell.scroll[0] = shell.scroll[1] = shell.scroll[2] = 0
  }
}

//Render stuff
function render(shell) {

  //Request next frame
  shell._rafHandle = requestAnimationFrame(shell._render)

  //Tick the shell
  tick(shell)
  
  //Compute frame time
  var dt
  if(shell._paused) {
    dt = shell._frameTime
  } else {
    dt = min(1.0, (hrtime() - shell._lastTick) / shell._tickRate)
  }
  
  //Draw a frame
  ++shell.frameCount
  var s = hrtime()
  shell.emit("render", dt)
  var t = hrtime()
  shell.frameTime = t - s
  
}

function isFocused(shell) {
  return (document.activeElement === document.body) ||
         (document.activeElement === shell.element)
}

//Set key up
function handleKeyUp(shell, ev) {
  ev.preventDefault()
  var kc = physicalKeyCode(ev.keyCode || ev.char || ev.which || ev.charCode)
  if(kc >= 0) {
    setKeyState(shell, kc, false)
  }
}

//Set key down
function handleKeyDown(shell, ev) {
  if(!isFocused(shell)) {
    return
  }
  if(ev.metaKey) {
    //Hack: Clear key state when meta gets pressed to prevent keys sticking
    handleBlur(shell, ev)
  } else {
    ev.preventDefault()
    var kc = physicalKeyCode(ev.keyCode || ev.char || ev.which || ev.charCode)
    if(kc >= 0) {
      setKeyState(shell, kc, true)
    }
  }
}

//Mouse events are really annoying
var mouseCodes = iota(32).map(function(n) {
  return virtualKeyCode("mouse-" + (n+1))
})

function setMouseButtons(shell, buttons) {
  for(var i=0; i<32; ++i) {
    setKeyState(shell, mouseCodes[i], !!(buttons & (1<<i)))
  }
}

function handleMouseMove(shell, ev) {
  if(shell._pointerLockActive) {
    var movementX = ev.movementX       ||
                    ev.mozMovementX    ||
                    ev.webkitMovementX ||
                    0,
        movementY = ev.movementY       ||
                    ev.mozMovementY    ||
                    ev.webkitMovementY ||
                    0
    shell.mouseX += movementX
    shell.mouseY += movementY
  } else {
    shell.mouseX = ev.clientX - shell.element.offsetLeft
    shell.mouseY = ev.clientY - shell.element.offsetTop
  }
  return false
}

function handleMouseDown(shell, ev) {
  setKeyState(shell, mouseCodes[ev.button], true)
  return false
}

function handleMouseUp(shell, ev) {
  setKeyState(shell, mouseCodes[ev.button], false)
  return false
}

function handleMouseEnter(shell, ev) {
  if(shell._pointerLockActive) {
    shell.prevMouseX = shell.mouseX = shell.width>>1
    shell.prevMouseY = shell.mouseY = shell.height>>1
  } else {
    shell.prevMouseX = shell.mouseX = ev.clientX - shell.element.offsetLeft
    shell.prevMouseY = shell.mouseY = ev.clientY - shell.element.offsetTop
  }
  return false
}

function handleMouseLeave(shell, ev) {
  setMouseButtons(shell, 0)
  return false
}

//Handle mouse wheel events
function handleMouseWheel(shell, ev) {
  var scale = 1
  switch(ev.deltaMode) {
    case 0: //Pixel
      scale = 1
    break
    case 1: //Line
      scale = 12
    break
    case 2: //Page
       scale = shell.height
    break
  }
  //Add scroll
  shell.scroll[0] +=  ev.deltaX * scale
  shell.scroll[1] +=  ev.deltaY * scale
  shell.scroll[2] += (ev.deltaZ * scale)||0.0
  return false
}

function handleContexMenu(shell, ev) {
  return false
}

function handleBlur(shell, ev) {
  var n = keyNames.length
    , c = shell._curKeyState
    , r = shell._releaseCount
    , i
  for(i=0; i<n; ++i) {
    if(c[i]) {
      ++r[i]
    }
    c[i] = false
  }
  return false
}

function handleResizeElement(shell, ev) {
  var w = shell.element.clientWidth|0
  var h = shell.element.clientHeight|0
  if((w !== shell._width) || (h !== shell._height)) {
    shell._width = w
    shell._height = h
    shell.emit("resize", w, h)
  }
}

function makeDefaultContainer() {
  var container = document.createElement("div")
  container.tabindex = 1
  container.style.position = "absolute"
  container.style.left = "0px"
  container.style.right = "0px"
  container.style.top = "0px"
  container.style.bottom = "0px"
  container.style.height = "100%"
  container.style.overflow = "hidden"
  document.body.appendChild(container)
  document.body.style.overflow = "hidden" //Prevent bounce
  document.body.style.height = "100%"
  return container
}

function createShell(options) {
  options = options || {}
  
  //Check fullscreen and pointer lock flags
  var useFullscreen = !!options.fullscreen
  var usePointerLock = useFullscreen
  if(typeof options.pointerLock !== undefined) {
    usePointerLock = !!options.pointerLock
  }
  
  //Create initial shell
  var shell = new GameShell()
  shell._tickRate = options.tickRate || 30
  shell.frameSkip = options.frameSkip || (shell._tickRate+5) * 5
  shell.stickyFullscreen = !!options.stickyFullscreen || !!options.sticky
  shell.stickyPointerLock = !!options.stickPointerLock || !options.sticky
  
  //Set bindings
  if(options.bindings) {
    shell.bindings = bindings
  }
  
  //Wait for dom to intiailize
  setTimeout(function() { domready(function initGameShell() {
    
    //Retrieve element
    var element = options.element
    if(typeof element === "string") {
      var e = document.querySelector(element)
      if(!e) {
        e = document.getElementById(element)
      }
      if(!e) {
        e = document.getElementByClass(element)[0]
      }
      if(!e) {
        e = makeDefaultContainer()
      }
      shell.element = e
    } else if(typeof element === "object" && !!element) {
      shell.element = element
    } else if(typeof element === "function") {
      shell.element = element()
    } else {
      shell.element = makeDefaultContainer()
    }
    
    //Disable user-select
    if(shell.element.style) {
      shell.element.style["-webkit-touch-callout"] = "none"
      shell.element.style["-webkit-user-select"] = "none"
      shell.element.style["-khtml-user-select"] = "none"
      shell.element.style["-moz-user-select"] = "none"
      shell.element.style["-ms-user-select"] = "none"
      shell.element.style["user-select"] = "none"
    }
    
    //Hook resize handler
    shell._width = shell.element.clientWidth
    shell._height = shell.element.clientHeight
    var handleResize = handleResizeElement.bind(undefined, shell)
    if(typeof MutationObserver !== "undefined") {
      var observer = new MutationObserver(handleResize)
      observer.observe(shell.element, {
        attributes: true,
        subtree: true
      })
    } else {
      shell.element.addEventListener("DOMSubtreeModified", handleResize, false)
    }
    window.addEventListener("resize", handleResize, false)
    
    //Hook keyboard listener
    window.addEventListener("keydown", handleKeyDown.bind(undefined, shell), false)
    window.addEventListener("keyup", handleKeyUp.bind(undefined, shell), false)
    
    //Disable right click
    shell.element.oncontextmenu = handleContexMenu.bind(undefined, shell)
    
    //Hook mouse listeners
    shell.element.addEventListener("mousedown", handleMouseDown.bind(undefined, shell), false)
    shell.element.addEventListener("mouseup", handleMouseUp.bind(undefined, shell), false)
    shell.element.addEventListener("mousemove", handleMouseMove.bind(undefined, shell), false)
    shell.element.addEventListener("mouseenter", handleMouseEnter.bind(undefined, shell), false)
    
    //Mouse leave
    var leave = handleMouseLeave.bind(undefined, shell)
    shell.element.addEventListener("mouseleave", leave, false)
    shell.element.addEventListener("mouseout", leave, false)
    window.addEventListener("mouseleave", leave, false)
    window.addEventListener("mouseout", leave, false)
    
    //Blur event 
    var blur = handleBlur.bind(undefined, shell)
    shell.element.addEventListener("blur", blur, false)
    shell.element.addEventListener("focusout", blur, false)
    shell.element.addEventListener("focus", blur, false)
    window.addEventListener("blur", blur, false)
    window.addEventListener("focusout", blur, false)
    window.addEventListener("focus", blur, false)

    //Mouse wheel handler
    addMouseWheel(shell.element, handleMouseWheel.bind(undefined, shell), false)

    //Fullscreen handler
    var fullscreenChange = handleFullscreen.bind(undefined, shell)
    document.addEventListener("fullscreenchange", fullscreenChange, false)
    document.addEventListener("mozfullscreenchange", fullscreenChange, false)
    document.addEventListener("webkitfullscreenchange", fullscreenChange, false)

    //Stupid fullscreen hack
    shell.element.addEventListener("click", tryFullscreen.bind(undefined, shell), false)

    //Pointer lock change handler
    var pointerLockChange = handlePointerLockChange.bind(undefined, shell)
    document.addEventListener("pointerlockchange", pointerLockChange, false)
    document.addEventListener("mozpointerlockchange", pointerLockChange, false)
    document.addEventListener("webkitpointerlockchange", pointerLockChange, false)
    document.addEventListener("pointerlocklost", pointerLockChange, false)
    document.addEventListener("webkitpointerlocklost", pointerLockChange, false)
    document.addEventListener("mozpointerlocklost", pointerLockChange, false)
    
    //Update flags
    shell.fullscreen = useFullscreen
    shell.pointerLock = usePointerLock
  
    //Default mouse button aliases
    shell.bind("mouse-left",   "mouse-1")
    shell.bind("mouse-right",  "mouse-3")
    shell.bind("mouse-middle", "mouse-2")
    
    //Initialize tick counter
    shell._lastTick = hrtime()
    shell.startTime = hrtime()

    //Unpause shell
    shell.paused = false
    
    //Emit initialize event
    shell.emit("init")
  })}, 0)
  
  return shell
}

module.exports = createShell

},{"./lib/hrtime-polyfill.js":"/Users/Romain/Documents/Programmation/Projets_Web/AntColony/node_modules/game-shell/lib/hrtime-polyfill.js","./lib/mousewheel-polyfill.js":"/Users/Romain/Documents/Programmation/Projets_Web/AntColony/node_modules/game-shell/lib/mousewheel-polyfill.js","./lib/raf-polyfill.js":"/Users/Romain/Documents/Programmation/Projets_Web/AntColony/node_modules/game-shell/lib/raf-polyfill.js","binary-search-bounds":"/Users/Romain/Documents/Programmation/Projets_Web/AntColony/node_modules/game-shell/node_modules/binary-search-bounds/search-bounds.js","domready":"/Users/Romain/Documents/Programmation/Projets_Web/AntColony/node_modules/game-shell/node_modules/domready/ready.js","events":"/usr/local/lib/node_modules/watchify/node_modules/browserify/node_modules/events/events.js","invert-hash":"/Users/Romain/Documents/Programmation/Projets_Web/AntColony/node_modules/game-shell/node_modules/invert-hash/invert.js","iota-array":"/Users/Romain/Documents/Programmation/Projets_Web/AntColony/node_modules/game-shell/node_modules/iota-array/iota.js","uniq":"/Users/Romain/Documents/Programmation/Projets_Web/AntColony/node_modules/game-shell/node_modules/uniq/uniq.js","util":"/usr/local/lib/node_modules/watchify/node_modules/browserify/node_modules/util/util.js","vkey":"/Users/Romain/Documents/Programmation/Projets_Web/AntColony/node_modules/game-shell/node_modules/vkey/index.js"}],"/Users/Romain/Documents/Programmation/Projets_Web/AntColony/node_modules/harmony-collections/harmony-collections.js":[function(require,module,exports){
/* (The MIT License)
 *
 * Copyright (c) 2012 Brandon Benvie <http://bbenvie.com>
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy of this software and
 * associated documentation files (the 'Software'), to deal in the Software without restriction,
 * including without limitation the rights to use, copy, modify, merge, publish, distribute,
 * sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included with all copies or
 * substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED 'AS IS', WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING
 * BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
 * NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY  CLAIM,
 * DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
 */

// Original WeakMap implementation by Gozala @ https://gist.github.com/1269991
// Updated and bugfixed by Raynos @ https://gist.github.com/1638059
// Expanded by Benvie @ https://github.com/Benvie/harmony-collections

void function(string_, object_, function_, prototype_, toString_,
              Array, Object, Function, FP, global, exports, undefined_, undefined){

  var getProperties = Object.getOwnPropertyNames,
      es5 = typeof getProperties === function_ && !(prototype_ in getProperties);

  var callbind = FP.bind
    ? FP.bind.bind(FP.call)
    : (function(call){
        return function(func){
          return function(){
            return call.apply(func, arguments);
          };
        };
      }(FP.call));

  var functionToString = callbind(FP[toString_]),
      objectToString = callbind({}[toString_]),
      numberToString = callbind(.0.toString),
      call = callbind(FP.call),
      apply = callbind(FP.apply),
      hasOwn = callbind({}.hasOwnProperty),
      push = callbind([].push),
      splice = callbind([].splice);

  var name = function(func){
    if (typeof func !== function_)
      return '';
    else if ('name' in func)
      return func.name;

    return functionToString(func).match(/^\n?function\s?(\w*)?_?\(/)[1];
  };

  var create = es5
    ? Object.create
    : function(proto, descs){
        var Ctor = function(){};
        Ctor[prototype_] = Object(proto);
        var object = new Ctor;

        if (descs)
          for (var key in descs)
            defineProperty(object, key, descs[k]);

        return object;
      };


  function Hash(){}

  if (es5) {
    void function(ObjectCreate){
      Hash.prototype = ObjectCreate(null);
      function inherit(obj){
        return ObjectCreate(obj);
      }
      Hash.inherit = inherit;
    }(Object.create);
  } else {
    void function(F){
      var iframe = document.createElement('iframe');
      iframe.style.display = 'none';
      document.body.appendChild(iframe);
      iframe.src = 'javascript:'
      Hash.prototype = iframe.contentWindow.Object.prototype;
      document.body.removeChild(iframe);
      iframe = null;

      var props = ['constructor', 'hasOwnProperty', 'propertyIsEnumerable',
                   'isProtoypeOf', 'toLocaleString', 'toString', 'valueOf'];

      for (var i=0; i < props.length; i++)
        delete Hash.prototype[props[i]];

      function inherit(obj){
        F.prototype = obj;
        obj = new F;
        F.prototype = null;
        return obj;
      }

      Hash.inherit = inherit;
    }(function(){});
  }

  var defineProperty = es5
    ? Object.defineProperty
    : function(object, key, desc) {
        object[key] = desc.value;
        return object;
      };

  var define = function(object, key, value){
    if (typeof key === function_) {
      value = key;
      key = name(value).replace(/_$/, '');
    }

    return defineProperty(object, key, { configurable: true, writable: true, value: value });
  };

  var isArray = es5
    ? (function(isArray){
        return function(o){
          return isArray(o) || o instanceof Array;
        };
      })(Array.isArray)
    : function(o){
        return o instanceof Array || objectToString(o) === '[object Array]';
      };

  // ############
  // ### Data ###
  // ############

  var builtinWeakMap = 'WeakMap' in global;

  var MapData = builtinWeakMap
    ? (function(){
      var BuiltinWeakMap = global.WeakMap,
          wmget = callbind(BuiltinWeakMap[prototype_].get),
          wmset = callbind(BuiltinWeakMap[prototype_].set),
          wmhas = callbind(BuiltinWeakMap[prototype_].has);

      function MapData(name){
        var map = new BuiltinWeakMap;

        this.get = function(o){
          return wmget(map, o);
        };
        this.set = function(o, v){
          wmset(map, o, v);
        };

        if (name) {
          this.wrap = function(o, v){
            if (wmhas(map, o))
              throw new TypeError("Object is already a " + name);
            wmset(map, o, v);
          };
          this.unwrap = function(o){
            var storage = wmget(map, o);
            if (!storage)
              throw new TypeError(name + " is not generic");
            return storage;
          };
        }
      }

      return MapData;
    })()
    : (function(){
      var locker = 'return function(k){if(k===s)return l}',
          random = Math.random,
          uids = new Hash,
          slice = callbind(''.slice),
          indexOf = callbind([].indexOf);

      var createUID = function(){
        var key = slice(numberToString(random(), 36), 2);
        return key in uids ? createUID() : uids[key] = key;
      };

      var globalID = createUID();

      // common per-object storage area made visible by patching getOwnPropertyNames'
      function getOwnPropertyNames(obj){
        var props = getProperties(obj);
        if (hasOwn(obj, globalID))
          splice(props, indexOf(props, globalID), 1);
        return props;
      }

      if (es5) {
        // check for the random key on an object, create new storage if missing, return it
        var storage = function(obj){
          if (!hasOwn(obj, globalID))
            defineProperty(obj, globalID, { value: new Hash });
          return obj[globalID];
        };

        define(Object, getOwnPropertyNames);
      } else {

        var toStringToString = function(s){
          function toString(){ return s }
          return toString[toString_] = toString;
        }(Object[prototype_][toString_]+'');

        // store the values on a custom valueOf in order to hide them but store them locally
        var storage = function(obj){
          if (hasOwn(obj, toString_) && globalID in obj[toString_])
            return obj[toString_][globalID];

          if (!(toString_ in obj))
            throw new Error("Can't store values for "+obj);

          var oldToString = obj[toString_];
          function toString(){ return oldToString.call(this) }
          obj[toString_] = toString;
          toString[toString_] = toStringToString;
          return toString[globalID] = {};
        };
      }



      // shim for [[MapData]] from es6 spec, and pulls double duty as WeakMap storage
      function MapData(name){
        var puid = createUID(),
            iuid = createUID(),
            secret = { value: undefined };

        var attach = function(obj){
          var store = storage(obj);
          if (hasOwn(store, puid))
            return store[puid](secret);

          var lockbox = new Hash;
          defineProperty(lockbox, iuid, secret);
          defineProperty(store, puid, {
            value: new Function('s', 'l', locker)(secret, lockbox)
          });
          return lockbox;
        };

        this.get = function(o){
          return attach(o)[iuid];
        };
        this.set = function(o, v){
          attach(o)[iuid] = v;
        };

        if (name) {
          this.wrap = function(o, v){
            var lockbox = attach(o);
            if (lockbox[iuid])
              throw new TypeError("Object is already a " + name);
            lockbox[iuid] = v;
          };
          this.unwrap = function(o){
            var storage = attach(o)[iuid];
            if (!storage)
              throw new TypeError(name + " is not generic");
            return storage;
          };
        }
      }

      return MapData;
    }());

  var exporter = (function(){
    // [native code] looks slightly different in each engine
    var src = (''+Object).split('Object');

    // fake [native code]
    function toString(){
      return src[0] + name(this) + src[1];
    }

    define(toString, toString);

    // attempt to use __proto__ so the methods don't all have an own toString
    var prepFunction = { __proto__: [] } instanceof Array
      ? function(func){ func.__proto__ = toString }
      : function(func){ define(func, toString) };

    // assemble an array of functions into a fully formed class
    var prepare = function(methods){
      var Ctor = methods.shift(),
          brand = '[object ' + name(Ctor) + ']';

      function toString(){ return brand }
      methods.push(toString);
      prepFunction(Ctor);

      for (var i=0; i < methods.length; i++) {
        prepFunction(methods[i]);
        define(Ctor[prototype_], methods[i]);
      }

      return Ctor;
    };

    return function(name, init){
      if (name in exports)
        return exports[name];

      var data = new MapData(name);

      return exports[name] = prepare(init(
        function(collection, value){
          data.wrap(collection, value);
        },
        function(collection){
          return data.unwrap(collection);
        }
      ));
    };
  }());


  // initialize collection with an iterable, currently only supports forEach function
  var initialize = function(iterable, callback){
    if (iterable !== null && typeof iterable === object_ && typeof iterable.forEach === function_) {
      iterable.forEach(function(item, i){
        if (isArray(item) && item.length === 2)
          callback(iterable[i][0], iterable[i][1]);
        else
          callback(iterable[i], i);
      });
    }
  }

  // attempt to fix the name of "delete_" methods, should work in V8 and spidermonkey
  var fixDelete = function(func, scopeNames, scopeValues){
    try {
      scopeNames[scopeNames.length] = ('return '+func).replace('e_', '\\u0065');
      return Function.apply(0, scopeNames).apply(0, scopeValues);
    } catch (e) {
      return func;
    }
  }

  var WM, HM, M;

  // ###############
  // ### WeakMap ###
  // ###############

  WM = builtinWeakMap ? (exports.WeakMap = global.WeakMap) : exporter('WeakMap', function(wrap, unwrap){
    var prototype = WeakMap[prototype_];
    var validate = function(key){
      if (key == null || typeof key !== object_ && typeof key !== function_)
        throw new TypeError("Invalid WeakMap key");
    };

    /**
     * @class        WeakMap
     * @description  Collection using objects with unique identities as keys that disallows enumeration
     *               and allows for better garbage collection.
     * @param        {Iterable} [iterable]  An item to populate the collection with.
     */
    function WeakMap(iterable){
      if (this === global || this == null || this === prototype)
        return new WeakMap(iterable);

      wrap(this, new MapData);

      var self = this;
      iterable && initialize(iterable, function(value, key){
        call(set, self, value, key);
      });
    }
    /**
     * @method       <get>
     * @description  Retrieve the value in the collection that matches key
     * @param        {Any} key
     * @return       {Any}
     */
    function get(key){
      validate(key);
      var value = unwrap(this).get(key);
      return value === undefined_ ? undefined : value;
    }
    /**
     * @method       <set>
     * @description  Add or update a pair in the collection. Enforces uniqueness by overwriting.
     * @param        {Any} key
     * @param        {Any} val
     **/
    function set(key, value){
      validate(key);
      // store a token for explicit undefined so that "has" works correctly
      unwrap(this).set(key, value === undefined ? undefined_ : value);
    }
    /*
     * @method       <has>
     * @description  Check if key is in the collection
     * @param        {Any} key
     * @return       {Boolean}
     **/
    function has(key){
      validate(key);
      return unwrap(this).get(key) !== undefined;
    }
    /**
     * @method       <delete>
     * @description  Remove key and matching value if found
     * @param        {Any} key
     * @return       {Boolean} true if item was in collection
     */
    function delete_(key){
      validate(key);
      var data = unwrap(this);

      if (data.get(key) === undefined)
        return false;

      data.set(key, undefined);
      return true;
    }

    delete_ = fixDelete(delete_, ['validate', 'unwrap'], [validate, unwrap]);
    return [WeakMap, get, set, has, delete_];
  });


  // ###############
  // ### HashMap ###
  // ###############

  HM = exporter('HashMap', function(wrap, unwrap){
    // separate numbers, strings, and atoms to compensate for key coercion to string

    var prototype = HashMap[prototype_],
        STRING = 0, NUMBER = 1, OTHER = 2,
        others = { 'true': true, 'false': false, 'null': null, 0: -0 };

    var proto = Math.random().toString(36).slice(2);

    var coerce = function(key){
      return key === '__proto__' ? proto : key;
    };

    var uncoerce = function(type, key){
      switch (type) {
        case STRING: return key === proto ? '__proto__' : key;
        case NUMBER: return +key;
        case OTHER: return others[key];
      }
    }


    var validate = function(key){
      if (key == null) return OTHER;
      switch (typeof key) {
        case 'boolean': return OTHER;
        case string_: return STRING;
        // negative zero has to be explicitly accounted for
        case 'number': return key === 0 && Infinity / key === -Infinity ? OTHER : NUMBER;
        default: throw new TypeError("Invalid HashMap key");
      }
    }

    /**
     * @class          HashMap
     * @description    Collection that only allows primitives to be keys.
     * @param          {Iterable} [iterable]  An item to populate the collection with.
     */
    function HashMap(iterable){
      if (this === global || this == null || this === prototype)
        return new HashMap(iterable);

      wrap(this, {
        size: 0,
        0: new Hash,
        1: new Hash,
        2: new Hash
      });

      var self = this;
      iterable && initialize(iterable, function(value, key){
        call(set, self, value, key);
      });
    }
    /**
     * @method       <get>
     * @description  Retrieve the value in the collection that matches key
     * @param        {Any} key
     * @return       {Any}
     */
    function get(key){
      return unwrap(this)[validate(key)][coerce(key)];
    }
    /**
     * @method       <set>
     * @description  Add or update a pair in the collection. Enforces uniqueness by overwriting.
     * @param        {Any} key
     * @param        {Any} val
     **/
    function set(key, value){
      var items = unwrap(this),
          data = items[validate(key)];

      key = coerce(key);
      key in data || items.size++;
      data[key] = value;
    }
    /**
     * @method       <has>
     * @description  Check if key exists in the collection.
     * @param        {Any} key
     * @return       {Boolean} is in collection
     **/
    function has(key){
      return coerce(key) in unwrap(this)[validate(key)];
    }
    /**
     * @method       <delete>
     * @description  Remove key and matching value if found
     * @param        {Any} key
     * @return       {Boolean} true if item was in collection
     */
    function delete_(key){
      var items = unwrap(this),
          data = items[validate(key)];

      key = coerce(key);
      if (key in data) {
        delete data[key];
        items.size--;
        return true;
      }

      return false;
    }
    /**
     * @method       <size>
     * @description  Retrieve the amount of items in the collection
     * @return       {Number}
     */
    function size(){
      return unwrap(this).size;
    }
    /**
     * @method       <forEach>
     * @description  Loop through the collection raising callback for each
     * @param        {Function} callback  `callback(value, key)`
     * @param        {Object}   context    The `this` binding for callbacks, default null
     */
    function forEach(callback, context){
      var data = unwrap(this);
      context = context == null ? global : context;
      for (var i=0; i < 3; i++)
        for (var key in data[i])
          call(callback, context, data[i][key], uncoerce(i, key), this);
    }

    delete_ = fixDelete(delete_, ['validate', 'unwrap', 'coerce'], [validate, unwrap, coerce]);
    return [HashMap, get, set, has, delete_, size, forEach];
  });


  // ###########
  // ### Map ###
  // ###########

  // if a fully implemented Map exists then use it
  if ('Map' in global && 'forEach' in global.Map.prototype) {
    M = exports.Map = global.Map;
  } else {
    M = exporter('Map', function(wrap, unwrap){
      // attempt to use an existing partially implemented Map
      var BuiltinMap = global.Map,
          prototype = Map[prototype_],
          wm = WM[prototype_],
          hm = (BuiltinMap || HM)[prototype_],
          mget    = [callbind(hm.get), callbind(wm.get)],
          mset    = [callbind(hm.set), callbind(wm.set)],
          mhas    = [callbind(hm.has), callbind(wm.has)],
          mdelete = [callbind(hm['delete']), callbind(wm['delete'])];

      var type = BuiltinMap
        ? function(){ return 0 }
        : function(o){ return +(typeof o === object_ ? o !== null : typeof o === function_) }

      // if we have a builtin Map we can let it do most of the heavy lifting
      var init = BuiltinMap
        ? function(){ return { 0: new BuiltinMap } }
        : function(){ return { 0: new HM, 1: new WM } };

      /**
       * @class         Map
       * @description   Collection that allows any kind of value to be a key.
       * @param         {Iterable} [iterable]  An item to populate the collection with.
       */
      function Map(iterable){
        if (this === global || this == null || this === prototype)
          return new Map(iterable);

        var data = init();
        data.keys = [];
        data.values = [];
        wrap(this, data);

        var self = this;
        iterable && initialize(iterable, function(value, key){
          call(set, self, value, key);
        });
      }
      /**
       * @method       <get>
       * @description  Retrieve the value in the collection that matches key
       * @param        {Any} key
       * @return       {Any}
       */
      function get(key){
        var data = unwrap(this),
            t = type(key);
        return data.values[mget[t](data[t], key)];
      }
      /**
       * @method       <set>
       * @description  Add or update a pair in the collection. Enforces uniqueness by overwriting.
       * @param        {Any} key
       * @param        {Any} val
       **/
      function set(key, value){
        var data = unwrap(this),
            t = type(key),
            index = mget[t](data[t], key);

        if (index === undefined) {
          mset[t](data[t], key, data.keys.length);
          push(data.keys, key);
          push(data.values, value);
        } else {
          data.keys[index] = key;
          data.values[index] = value;
        }
      }
      /**
       * @method       <has>
       * @description  Check if key exists in the collection.
       * @param        {Any} key
       * @return       {Boolean} is in collection
       **/
      function has(key){
        var t = type(key);
        return mhas[t](unwrap(this)[t], key);
      }
      /**
       * @method       <delete>
       * @description  Remove key and matching value if found
       * @param        {Any} key
       * @return       {Boolean} true if item was in collection
       */
      function delete_(key){
        var data = unwrap(this),
            t = type(key),
            index = mget[t](data[t], key);

        if (index === undefined)
          return false;

        mdelete[t](data[t], key);
        splice(data.keys, index, 1);
        splice(data.values, index, 1);
        return true;
      }
      /**
       * @method       <size>
       * @description  Retrieve the amount of items in the collection
       * @return       {Number}
       */
      function size(){
        return unwrap(this).keys.length;
      }
      /**
       * @method       <forEach>
       * @description  Loop through the collection raising callback for each
       * @param        {Function} callback  `callback(value, key)`
       * @param        {Object}   context    The `this` binding for callbacks, default null
       */
      function forEach(callback, context){
        var data = unwrap(this),
            keys = data.keys,
            values = data.values;

        context = context == null ? global : context;

        for (var i=0, len=keys.length; i < len; i++)
          call(callback, context, values[i], keys[i], this);
      }

      delete_ = fixDelete(delete_,
        ['type', 'unwrap', 'call', 'splice'],
        [type, unwrap, call, splice]
      );
      return [Map, get, set, has, delete_, size, forEach];
    });
  }


  // ###########
  // ### Set ###
  // ###########

  exporter('Set', function(wrap, unwrap){
    var prototype = Set[prototype_],
        m = M[prototype_],
        msize = callbind(m.size),
        mforEach = callbind(m.forEach),
        mget = callbind(m.get),
        mset = callbind(m.set),
        mhas = callbind(m.has),
        mdelete = callbind(m['delete']);

    /**
     * @class        Set
     * @description  Collection of values that enforces uniqueness.
     * @param        {Iterable} [iterable]  An item to populate the collection with.
     **/
    function Set(iterable){
      if (this === global || this == null || this === prototype)
        return new Set(iterable);

      wrap(this, new M);

      var self = this;
      iterable && initialize(iterable, function(value, key){
        call(add, self, key);
      });
    }
    /**
     * @method       <add>
     * @description  Insert value if not found, enforcing uniqueness.
     * @param        {Any} val
     */
    function add(key){
      mset(unwrap(this), key, key);
    }
    /**
     * @method       <has>
     * @description  Check if key exists in the collection.
     * @param        {Any} key
     * @return       {Boolean} is in collection
     **/
    function has(key){
      return mhas(unwrap(this), key);
    }
    /**
     * @method       <delete>
     * @description  Remove key and matching value if found
     * @param        {Any} key
     * @return       {Boolean} true if item was in collection
     */
    function delete_(key){
      return mdelete(unwrap(this), key);
    }
    /**
     * @method       <size>
     * @description  Retrieve the amount of items in the collection
     * @return       {Number}
     */
    function size(){
      return msize(unwrap(this));
    }
    /**
     * @method       <forEach>
     * @description  Loop through the collection raising callback for each. Index is simply the counter for the current iteration.
     * @param        {Function} callback  `callback(value, index)`
     * @param        {Object}   context    The `this` binding for callbacks, default null
     */
    function forEach(callback, context){
      var index = 0,
          self = this;
      mforEach(unwrap(this, function(key){
        call(callback, this, key, index++, self);
      }, context));
    }

    delete_ = fixDelete(delete_, ['mdelete', 'unwrap'], [mdelete, unwrap]);
    return [Set, add, has, delete_, size, forEach];
  });
}('string', 'object', 'function', 'prototype', 'toString',
  Array, Object, Function, Function.prototype, (0, eval)('this'),
  typeof exports === 'undefined' ? this : exports, {});

},{}],"/Users/Romain/Documents/Programmation/Projets_Web/AntColony/node_modules/parse-svg-path/index.js":[function(require,module,exports){

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

},{}],"/Users/Romain/Documents/Programmation/Projets_Web/AntColony/src/index.js":[function(require,module,exports){
"use strict"

var shell = require("game-shell")();
var Map = require("harmony-collections").Map;
var Set = require("harmony-collections").Set;
var dt = require("delaunay-triangulate");
var parse = require('parse-svg-path')


var nbRandomPoints = 100;
var nbAnts = 200;
var textMesh = false;
var nbCity = 10;

var sqrt = Math.sqrt;
var pow = Math.pow;
var floor = Math.floor;
var random = Math.random;

var svgString = "m 1246.3864,422.67046 c -32.6566,-3.05977 -80.9202,-1.60892 -64,-49.23572 -1.4759,-17.26999 -2.9389,-34.52697 19.6179,-27.76428 28.2051,-7.71782 5.6299,36.53428 29.3821,44 15.6401,11.08789 36.8268,15.16788 55.7601,11.43591 18.9333,-3.73198 35.6132,-15.27593 42.2399,-35.43591 11.114,-30.53447 -5.455,-63.66997 -34.7792,-76 -21.2379,-12.18874 -45.1963,-24.31647 -65.5186,-39.83315 -20.3223,-15.51668 -37.0085,-34.42231 -43.7022,-60.16685 -5.2777,-23.82651 1.2957,-47.48151 15.0269,-66.10238 13.7312,-18.62088 34.6201,-32.207623 57.9731,-35.897619 19.5714,-4.408067 39.5199,-4.455217 59.3132,-2.055584 19.7933,2.399634 39.4313,7.246051 58.3817,12.62512 1.5096,19.173023 5.7329,49.670133 -0.3128,65.430463 -16.657,-0.26915 -39.1842,5.88414 -31.6548,-19 -1.6482,-23.82053 -21.3427,-36.82701 -42.4938,-38.50012 -21.1511,-1.67312 -43.7588,7.98715 -51.2335,29.50012 -14.1826,27.84603 4.1914,57.9202 30,70 19.9649,13.49423 43.5526,24.76438 64.3552,38.77589 20.8027,14.01152 38.8202,30.7644 47.6448,55.22411 8.8577,24.73667 4.7321,50.90082 -7.7951,72.42372 -12.5272,21.5229 -33.456,38.40456 -58.2049,44.57628 -25.5516,8.81359 -53.4285,7.61391 -80,6 z m -502.75004,-4.49583 c -14.49634,-20.16798 -29.05177,-40.31989 -43.66494,-60.44416 -14.61316,-20.12428 -29.28404,-40.22091 -44.01126,-60.27834 -14.72723,-20.05744 -29.51079,-40.07567 -44.34931,-60.04313 -14.83851,-19.96747 -29.73199,-39.88417 -44.67904,-59.73854 0.085,67.94114 -2.59454,136.28635 1.45455,204 2.37361,21.33297 41.95968,2.60802 32,29.23572 -5.68484,10.73778 -30.97918,1.86608 -44.5036,4.76428 -19.16547,0 -38.33093,0 -57.4964,0 -12.17745,-31.99983 32.34126,-11.09351 34,-38 0.9227,-19.95147 1.50618,-39.93087 1.84747,-59.92634 0.34129,-19.99547 0.44039,-40.00702 0.39435,-60.02281 -0.0461,-20.01579 -0.23725,-40.03582 -0.47657,-60.04824 -0.23931,-20.01242 -0.52674,-40.01724 -0.76525,-60.00261 10.20299,-32.81396 -63.88039,-18.66932 -25.06566,-45.731739 23.84724,1.313312 50.33669,-2.975416 73.00768,1.582846 14.56739,20.242453 29.23861,40.419643 43.95531,60.568633 14.71671,20.14899 29.47889,40.26977 44.2282,60.39943 14.74931,20.12965 29.48574,40.26817 44.15094,60.45262 14.66519,20.18445 29.25916,40.41483 43.72353,60.72821 0.83332,-70.07002 2.82145,-140.78692 -1,-210.54546 -10.39601,-9.71246 -51.54014,-21.85965 -24.14642,-33.454539 30.9003,0.49167 61.92696,-1.293841 92.71785,1.42857 12.44918,33.638319 -51.66003,12.223739 -34.00692,55.415179 -0.24121,22.59532 -0.38824,45.19117 -0.47532,67.78734 -0.0871,22.59618 -0.11419,45.19267 -0.11559,67.7893 -10e-4,22.59662 0.0229,45.19336 0.0387,67.79003 0.0158,22.59667 0.0231,45.19326 -0.0123,67.78958 -15.5356,-0.53883 -31.43722,1.23098 -46.75,-1.49583 z M 132.38637,405.67046 c -2.30195,-19.38098 35.20853,-6.01885 35.99999,-32 10.6943,-23.20341 21.35806,-46.41833 32.00899,-69.63818 10.65093,-23.21985 21.28905,-46.44461 31.93207,-69.6677 10.64303,-23.22308 21.29096,-46.44449 31.96154,-69.65761 10.67058,-23.21313 21.3638,-46.41797 32.0974,-69.607939 14.00633,-14.784664 31.06301,-3.589508 34,14.571429 9.47589,20.50998 18.90267,41.04032 28.32168,61.57386 9.41901,20.53355 18.83027,41.07031 28.27512,61.59315 9.44485,20.52284 18.92331,41.03177 28.47672,61.50964 9.55341,20.47788 19.18179,40.9247 28.92648,61.32335 8.94847,18.57494 17.97684,44.59976 44,40 6.95283,23.21831 -11.22619,21.30477 -28.79645,20 -17.53393,0 -35.06785,0 -52.60178,0 -17.53392,0 -35.06785,0 -52.60177,0 -13.61053,-31.56953 34.86478,-11.43099 30,-37 -7.96326,-20.44937 -16.90067,-40.88219 -27.17188,-60.1481 -22.69835,-0.30159 -45.74783,-1.00478 -68.70149,-1.13746 -22.95367,-0.13269 -45.81153,0.30513 -68.12663,2.28556 -10.15614,23.48394 -27.65652,46.62921 -26.57143,72.57143 21.87983,-9.59529 37.44237,23.82765 14.09999,23.42857 -25.17618,0 -50.35237,0 -75.52855,0 0,-3.33333 0,-6.66667 0,-10 z m 213.99999,-112 c -9.85189,-19.92913 -19.12968,-40.10864 -28.4649,-60.26037 -9.33521,-20.15173 -18.72785,-40.27567 -28.8094,-60.09366 -8.0471,4.62582 -13.12108,26.97607 -19.7257,38.35403 -12.60829,27.21944 -25.26661,54.41816 -37,82 18.19346,1.65846 37.60009,2.48757 57.00505,2.48751 19.40495,-6e-5 38.80822,-0.82929 56.99495,-2.48751 z m 582,112 c -1.2337,-21.88184 52.78642,-3.94103 39.43532,-42.15625 0.24195,-20.8192 0.38925,-41.63897 0.47635,-62.45909 0.0871,-20.82012 0.11398,-41.64058 0.1151,-62.46118 10e-4,-20.8206 -0.0235,-41.64134 -0.0394,-62.46199 -0.0159,-20.82065 -0.0232,-41.64122 0.0127,-62.46149 -26.37832,2.29828 -55.93232,-6.20371 -79.85787,6.80589 -1.50256,28.04473 -4.68344,46.64651 -35.23142,39.19411 -7.11643,-15.54116 -0.90104,-46.83584 -1.48214,-66.571429 25.02495,-0.70812 50.0617,-1.115779 75.10565,-1.337873 25.04395,-0.222093 50.0951,-0.258621 75.14885,-0.224476 25.0538,0.03414 50.1101,0.138961 75.1645,0.199555 25.0544,0.06059 50.1067,0.07697 75.1525,-0.06578 0,22.666673 0,45.333333 0,68.000003 -11.5708,-3.09584 -38.1128,8.77945 -32,-12 6.3763,-41.53883 -40.2943,-34.36261 -67.7072,-34 -23.6002,-7.69784 -14.9407,16.86409 -16.2929,31.41431 0.1745,19.75739 0.1467,39.53024 0.067,59.30763 -0.08,19.77738 -0.2116,39.5593 -0.2452,59.33481 -0.034,19.77551 0.031,39.54462 0.3439,59.29638 0.3129,19.75176 0.8743,39.48618 1.8343,59.19232 6.3204,22.76756 50.3983,0.0574 38.0001,33.45455 -24.6667,0 -49.3333,0 -74,0 -24.66671,0 -49.3334,0 -74.00007,0 0,-3.33333 0,-6.66667 0,-10 z";

function svgToPoints(svgString) {
    var points = [];
    var edges = [];

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
                break;
            case "M":
                X = command[1];
                Y = command[2];
                prevPoint = undefined;
                break;  
            case "c":
                X += command[5];
                Y += command[6];
                points.push({id:nbPoints, x:X, y:Y});
                nbPoints++;
                if (prevPoint) {
                    edges.push([prevPoint, nbPoints]);
                    prevPoint = nbPoints;
                }
                break;    
        }
    }
    return {points : points, edges : edges};
}

function range(start, count) {
    return Array.apply(0, Array(count)).map(function (element, index) { return index + start });
}

function sign(x) { return x ? x < 0 ? -1 : 1 : 0; }

// initialize points
var points = [];
var citySet;

if (textMesh){

    var myText = svgToPoints(svgString);
    points = myText.points;
    citySet = new Set(range(0, points.length));
    var scale = 0.5

    // scale points to [0,1] + scale
    var maxX = Math.max.apply(Math, points.map(function(p){return p.x}));
    var minX = Math.min.apply(Math, points.map(function(p){return p.x}));
    var maxY = Math.max.apply(Math, points.map(function(p){return p.y}));
    var minY = Math.min.apply(Math, points.map(function(p){return p.y}));
    points = points.map(function(p){return {id:p.id, x: 0.4*(p.x-minX)/(maxX-minX)+0.25, y: 0.4*(p.y-minY)/(maxY-minY)+0.25}})

    //add random points
    var nbPoints = points.length;
    for(var i=0; i<nbRandomPoints; ++i) {
        points.push({id : nbPoints, x:random(), y:random()});
        nbPoints++;
    }

} else {
    //add random points
    var nbPoints = 0;
    for(var i=0; i<nbRandomPoints; ++i) {
        points.push({id : nbPoints, x:random(), y:random()});
        nbPoints++;
    }
    citySet = new Set(range(0, nbCity));
}




// triangulate
var cells = dt(points.map(function(p){return [p.x, p.y]}))

// create edges
var nextEdges = new Map();
var edges = [];
var permutations = [[0,1], [1,0], [0,2], [2,0], [1,2], [2,1]];
var nbEdges = 0;
cells.forEach(function(cell){
  for(var i=0; i<permutations.length; ++i){
    var s = permutations[i][0];
    var d = permutations[i][1];
    var ps = points[cell[s]];
    var pd = points[cell[d]];
    var distance = sqrt( pow(ps.x - pd.x, 2) + pow(ps.y - pd.y, 2) );
    var edge = {id : nbEdges,
          source: cell[s], 
          destination: cell[d], 
          distance : distance,
          direction : Math.atan((pd.y-ps.y)/(pd.x-ps.x)),
          orientation : sign((pd.x-ps.x)),
          pheromon : 1/distance
          };
    var nexts;
    if(nextEdges.has(ps.id)){
        nexts = nextEdges.get(ps.id);
        nexts.push(edge);
    } else {
        nexts = [edge];
    }
    nextEdges.set(ps.id, nexts);
    edges.push(edge);
    nbEdges++;
  }
  
})

// initialize ants
var population = new Array(nbAnts);
var i,j;
for (i = 0; i < nbAnts; i++) {
    // take a random edge
    var edge = edges[Math.floor(edges.length*random())];
    var x = points[edge.source].x 
    var y = points[edge.source].y
    population[i] = new Ant(x, y, edge);
}

var canvas, context;

shell.on("init", function() {
    canvas = document.createElement("canvas");
    canvas.width = shell.width;
    canvas.height = shell.height;
    context = canvas.getContext("2d");
    shell.element.appendChild(canvas);
})

shell.on("resize", function(w, h) {
    canvas.width = w;
    canvas.height = h;
})

shell.on("render", function() {
    var w = canvas.width;
    var h = canvas.height;
    var mouse = [shell.mouseX/w, shell.mouseY/h];
    context.setTransform(w, 0, 0, h, 0, 0);
    context.fillStyle = "#fff";
    context.fillRect(0,0,w,h);

    // edges
    context.strokeStyle = "#000";
    for(var i=0; i<edges.length; ++i) {
        var edge = edges[i];
        if (edge.pheromon != 0){
            context.lineWidth = 0.00001 * edge.pheromon;
        }else {
            context.lineWidth = 0.00001;
        }
        context.beginPath();
        context.moveTo(points[edge.source].x, points[edge.source].y);
        context.lineTo(points[edge.destination].x, points[edge.destination].y);
        context.stroke();
    }

    // vertices
    for(var i=0; i<points.length; ++i) {
        context.beginPath()
        var point = points[i];
        if (citySet.has(point.id)) {
            context.fillStyle = "#0101DF";
            context.arc(point.x, point.y, 0.006, 0, 2*Math.PI);
        }
        else {
            context.fillStyle = "#000";
            context.arc(points[i].x, points[i].y, 0.003, 0, 2*Math.PI);
        }
        context.closePath();
        context.fill();
    }

    // move ants
    for (i = 0; i < nbAnts; i++) {
        population[i].transit();
    }

    // pheromon evaporation
    for (i = 0; i < edges.length; i++) {
        if(edges[i].pheromon > 0){
            edges[i].pheromon -= 0.001;
        }
    }


    for(var i=0; i<population.length; ++i) {
        context.beginPath()
        var x = population[i].posX //+ 0.01*random();
        var y = population[i].posY //+ 0.01*random();
        if (population[i].state === "pheromoning"){context.fillStyle = "#FF0000"}
        else {context.fillStyle = "#610B0B"}
        context.arc(x, y, 0.003, 0, 2*Math.PI)
        context.closePath()
        context.fill()
    }
  
})

function Ant(x, y, edge) {                                            
    this.posX = x;                
    this.posY = y;
    this.edge = edge;
    this.step = 0;
    this.state = "forage";
    this.transit = statemachine; 
    this.move = move;
    this.edges = [];
    this.lastCity = undefined;
}
// forage: the ant wanders around without any pheromon deposition
// once it finds a city, it starts remembering the nodes it goes through
// when it finds another city, it computes the path length and adds pheromons one each edges
// proportionnaly to the shortestness of the path
// it resets the list of nodes and continues
// while foraging the ant choses the path with a pheromon preference

function move() {
    var edgeChanged;
    var cityReached = false;
    // on edge
    if (this.step < this.edge.distance){
        this.posX += 0.005*Math.cos(this.edge.direction)*this.edge.orientation;
        this.posY += 0.005*Math.sin(this.edge.direction)*this.edge.orientation;
        this.step += 0.005;
        edgeChanged = false;
    // on vertex
    } else {
        this.step = 0;
        this.posX = points[this.edge.destination].x;
        this.posY = points[this.edge.destination].y;
        var possibleEdges = nextEdges.get(this.edge.destination);
        // flip a coin and either take the smelliest path of a random one
        if (random() > 0.5){
            var smells = possibleEdges.map(function(e){return e.pheromon});
            var index = smells.indexOf(Math.max.apply(Math, smells));
            this.edge = possibleEdges[index];
        } else {
            this.edge = possibleEdges[floor(random()*possibleEdges.length)];
        }
        cityReached = citySet.has(this.edge.source);
        edgeChanged = true;
    }
    return {cityReached: cityReached, edgeChanged: edgeChanged};
}


function statemachine() {
    switch (this.state) {
        case "forage":
            var res = this.move();
            if (res.cityReached) {
                this.state = "pheromoning";
                this.lastCity = this.edge.source;
            };
            break;
        case "pheromoning":
            var res = this.move();
            if (res.edgeChanged) {
                this.edges.push(this.edge);
                // found a city
                if (res.cityReached && (this.edge.source != this.lastCity) ){
                    // compute the length of the path
                    var pathLength = this.edges.map(function(e){return e.distance}).reduce(function(a,b){return a + b});
                    var deltaPheromone = 1/pathLength;
                    this.edges.forEach(function(e){e.pheromon += deltaPheromone});
                    // console.log(deltaPheromone, this.edges);
                    this.edges = [this.edge];
                    this.lastCity = this.edge.source;
                }
            }
          break;

    }
}
},{"delaunay-triangulate":"/Users/Romain/Documents/Programmation/Projets_Web/AntColony/node_modules/delaunay-triangulate/triangulate.js","game-shell":"/Users/Romain/Documents/Programmation/Projets_Web/AntColony/node_modules/game-shell/shell.js","harmony-collections":"/Users/Romain/Documents/Programmation/Projets_Web/AntColony/node_modules/harmony-collections/harmony-collections.js","parse-svg-path":"/Users/Romain/Documents/Programmation/Projets_Web/AntColony/node_modules/parse-svg-path/index.js"}],"/usr/local/lib/node_modules/watchify/node_modules/browserify/node_modules/events/events.js":[function(require,module,exports){
// Copyright Joyent, Inc. and other Node contributors.
//
// Permission is hereby granted, free of charge, to any person obtaining a
// copy of this software and associated documentation files (the
// "Software"), to deal in the Software without restriction, including
// without limitation the rights to use, copy, modify, merge, publish,
// distribute, sublicense, and/or sell copies of the Software, and to permit
// persons to whom the Software is furnished to do so, subject to the
// following conditions:
//
// The above copyright notice and this permission notice shall be included
// in all copies or substantial portions of the Software.
//
// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS
// OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
// MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN
// NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM,
// DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR
// OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE
// USE OR OTHER DEALINGS IN THE SOFTWARE.

function EventEmitter() {
  this._events = this._events || {};
  this._maxListeners = this._maxListeners || undefined;
}
module.exports = EventEmitter;

// Backwards-compat with node 0.10.x
EventEmitter.EventEmitter = EventEmitter;

EventEmitter.prototype._events = undefined;
EventEmitter.prototype._maxListeners = undefined;

// By default EventEmitters will print a warning if more than 10 listeners are
// added to it. This is a useful default which helps finding memory leaks.
EventEmitter.defaultMaxListeners = 10;

// Obviously not all Emitters should be limited to 10. This function allows
// that to be increased. Set to zero for unlimited.
EventEmitter.prototype.setMaxListeners = function(n) {
  if (!isNumber(n) || n < 0 || isNaN(n))
    throw TypeError('n must be a positive number');
  this._maxListeners = n;
  return this;
};

EventEmitter.prototype.emit = function(type) {
  var er, handler, len, args, i, listeners;

  if (!this._events)
    this._events = {};

  // If there is no 'error' event listener then throw.
  if (type === 'error') {
    if (!this._events.error ||
        (isObject(this._events.error) && !this._events.error.length)) {
      er = arguments[1];
      if (er instanceof Error) {
        throw er; // Unhandled 'error' event
      } else {
        throw TypeError('Uncaught, unspecified "error" event.');
      }
      return false;
    }
  }

  handler = this._events[type];

  if (isUndefined(handler))
    return false;

  if (isFunction(handler)) {
    switch (arguments.length) {
      // fast cases
      case 1:
        handler.call(this);
        break;
      case 2:
        handler.call(this, arguments[1]);
        break;
      case 3:
        handler.call(this, arguments[1], arguments[2]);
        break;
      // slower
      default:
        len = arguments.length;
        args = new Array(len - 1);
        for (i = 1; i < len; i++)
          args[i - 1] = arguments[i];
        handler.apply(this, args);
    }
  } else if (isObject(handler)) {
    len = arguments.length;
    args = new Array(len - 1);
    for (i = 1; i < len; i++)
      args[i - 1] = arguments[i];

    listeners = handler.slice();
    len = listeners.length;
    for (i = 0; i < len; i++)
      listeners[i].apply(this, args);
  }

  return true;
};

EventEmitter.prototype.addListener = function(type, listener) {
  var m;

  if (!isFunction(listener))
    throw TypeError('listener must be a function');

  if (!this._events)
    this._events = {};

  // To avoid recursion in the case that type === "newListener"! Before
  // adding it to the listeners, first emit "newListener".
  if (this._events.newListener)
    this.emit('newListener', type,
              isFunction(listener.listener) ?
              listener.listener : listener);

  if (!this._events[type])
    // Optimize the case of one listener. Don't need the extra array object.
    this._events[type] = listener;
  else if (isObject(this._events[type]))
    // If we've already got an array, just append.
    this._events[type].push(listener);
  else
    // Adding the second element, need to change to array.
    this._events[type] = [this._events[type], listener];

  // Check for listener leak
  if (isObject(this._events[type]) && !this._events[type].warned) {
    var m;
    if (!isUndefined(this._maxListeners)) {
      m = this._maxListeners;
    } else {
      m = EventEmitter.defaultMaxListeners;
    }

    if (m && m > 0 && this._events[type].length > m) {
      this._events[type].warned = true;
      console.error('(node) warning: possible EventEmitter memory ' +
                    'leak detected. %d listeners added. ' +
                    'Use emitter.setMaxListeners() to increase limit.',
                    this._events[type].length);
      if (typeof console.trace === 'function') {
        // not supported in IE 10
        console.trace();
      }
    }
  }

  return this;
};

EventEmitter.prototype.on = EventEmitter.prototype.addListener;

EventEmitter.prototype.once = function(type, listener) {
  if (!isFunction(listener))
    throw TypeError('listener must be a function');

  var fired = false;

  function g() {
    this.removeListener(type, g);

    if (!fired) {
      fired = true;
      listener.apply(this, arguments);
    }
  }

  g.listener = listener;
  this.on(type, g);

  return this;
};

// emits a 'removeListener' event iff the listener was removed
EventEmitter.prototype.removeListener = function(type, listener) {
  var list, position, length, i;

  if (!isFunction(listener))
    throw TypeError('listener must be a function');

  if (!this._events || !this._events[type])
    return this;

  list = this._events[type];
  length = list.length;
  position = -1;

  if (list === listener ||
      (isFunction(list.listener) && list.listener === listener)) {
    delete this._events[type];
    if (this._events.removeListener)
      this.emit('removeListener', type, listener);

  } else if (isObject(list)) {
    for (i = length; i-- > 0;) {
      if (list[i] === listener ||
          (list[i].listener && list[i].listener === listener)) {
        position = i;
        break;
      }
    }

    if (position < 0)
      return this;

    if (list.length === 1) {
      list.length = 0;
      delete this._events[type];
    } else {
      list.splice(position, 1);
    }

    if (this._events.removeListener)
      this.emit('removeListener', type, listener);
  }

  return this;
};

EventEmitter.prototype.removeAllListeners = function(type) {
  var key, listeners;

  if (!this._events)
    return this;

  // not listening for removeListener, no need to emit
  if (!this._events.removeListener) {
    if (arguments.length === 0)
      this._events = {};
    else if (this._events[type])
      delete this._events[type];
    return this;
  }

  // emit removeListener for all listeners on all events
  if (arguments.length === 0) {
    for (key in this._events) {
      if (key === 'removeListener') continue;
      this.removeAllListeners(key);
    }
    this.removeAllListeners('removeListener');
    this._events = {};
    return this;
  }

  listeners = this._events[type];

  if (isFunction(listeners)) {
    this.removeListener(type, listeners);
  } else {
    // LIFO order
    while (listeners.length)
      this.removeListener(type, listeners[listeners.length - 1]);
  }
  delete this._events[type];

  return this;
};

EventEmitter.prototype.listeners = function(type) {
  var ret;
  if (!this._events || !this._events[type])
    ret = [];
  else if (isFunction(this._events[type]))
    ret = [this._events[type]];
  else
    ret = this._events[type].slice();
  return ret;
};

EventEmitter.listenerCount = function(emitter, type) {
  var ret;
  if (!emitter._events || !emitter._events[type])
    ret = 0;
  else if (isFunction(emitter._events[type]))
    ret = 1;
  else
    ret = emitter._events[type].length;
  return ret;
};

function isFunction(arg) {
  return typeof arg === 'function';
}

function isNumber(arg) {
  return typeof arg === 'number';
}

function isObject(arg) {
  return typeof arg === 'object' && arg !== null;
}

function isUndefined(arg) {
  return arg === void 0;
}

},{}],"/usr/local/lib/node_modules/watchify/node_modules/browserify/node_modules/inherits/inherits_browser.js":[function(require,module,exports){
if (typeof Object.create === 'function') {
  // implementation from standard node.js 'util' module
  module.exports = function inherits(ctor, superCtor) {
    ctor.super_ = superCtor
    ctor.prototype = Object.create(superCtor.prototype, {
      constructor: {
        value: ctor,
        enumerable: false,
        writable: true,
        configurable: true
      }
    });
  };
} else {
  // old school shim for old browsers
  module.exports = function inherits(ctor, superCtor) {
    ctor.super_ = superCtor
    var TempCtor = function () {}
    TempCtor.prototype = superCtor.prototype
    ctor.prototype = new TempCtor()
    ctor.prototype.constructor = ctor
  }
}

},{}],"/usr/local/lib/node_modules/watchify/node_modules/browserify/node_modules/process/browser.js":[function(require,module,exports){
// shim for using process in browser

var process = module.exports = {};

process.nextTick = (function () {
    var canSetImmediate = typeof window !== 'undefined'
    && window.setImmediate;
    var canPost = typeof window !== 'undefined'
    && window.postMessage && window.addEventListener
    ;

    if (canSetImmediate) {
        return function (f) { return window.setImmediate(f) };
    }

    if (canPost) {
        var queue = [];
        window.addEventListener('message', function (ev) {
            var source = ev.source;
            if ((source === window || source === null) && ev.data === 'process-tick') {
                ev.stopPropagation();
                if (queue.length > 0) {
                    var fn = queue.shift();
                    fn();
                }
            }
        }, true);

        return function nextTick(fn) {
            queue.push(fn);
            window.postMessage('process-tick', '*');
        };
    }

    return function nextTick(fn) {
        setTimeout(fn, 0);
    };
})();

process.title = 'browser';
process.browser = true;
process.env = {};
process.argv = [];

function noop() {}

process.on = noop;
process.addListener = noop;
process.once = noop;
process.off = noop;
process.removeListener = noop;
process.removeAllListeners = noop;
process.emit = noop;

process.binding = function (name) {
    throw new Error('process.binding is not supported');
}

// TODO(shtylman)
process.cwd = function () { return '/' };
process.chdir = function (dir) {
    throw new Error('process.chdir is not supported');
};

},{}],"/usr/local/lib/node_modules/watchify/node_modules/browserify/node_modules/util/support/isBufferBrowser.js":[function(require,module,exports){
module.exports = function isBuffer(arg) {
  return arg && typeof arg === 'object'
    && typeof arg.copy === 'function'
    && typeof arg.fill === 'function'
    && typeof arg.readUInt8 === 'function';
}
},{}],"/usr/local/lib/node_modules/watchify/node_modules/browserify/node_modules/util/util.js":[function(require,module,exports){
(function (process,global){
// Copyright Joyent, Inc. and other Node contributors.
//
// Permission is hereby granted, free of charge, to any person obtaining a
// copy of this software and associated documentation files (the
// "Software"), to deal in the Software without restriction, including
// without limitation the rights to use, copy, modify, merge, publish,
// distribute, sublicense, and/or sell copies of the Software, and to permit
// persons to whom the Software is furnished to do so, subject to the
// following conditions:
//
// The above copyright notice and this permission notice shall be included
// in all copies or substantial portions of the Software.
//
// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS
// OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
// MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN
// NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM,
// DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR
// OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE
// USE OR OTHER DEALINGS IN THE SOFTWARE.

var formatRegExp = /%[sdj%]/g;
exports.format = function(f) {
  if (!isString(f)) {
    var objects = [];
    for (var i = 0; i < arguments.length; i++) {
      objects.push(inspect(arguments[i]));
    }
    return objects.join(' ');
  }

  var i = 1;
  var args = arguments;
  var len = args.length;
  var str = String(f).replace(formatRegExp, function(x) {
    if (x === '%%') return '%';
    if (i >= len) return x;
    switch (x) {
      case '%s': return String(args[i++]);
      case '%d': return Number(args[i++]);
      case '%j':
        try {
          return JSON.stringify(args[i++]);
        } catch (_) {
          return '[Circular]';
        }
      default:
        return x;
    }
  });
  for (var x = args[i]; i < len; x = args[++i]) {
    if (isNull(x) || !isObject(x)) {
      str += ' ' + x;
    } else {
      str += ' ' + inspect(x);
    }
  }
  return str;
};


// Mark that a method should not be used.
// Returns a modified function which warns once by default.
// If --no-deprecation is set, then it is a no-op.
exports.deprecate = function(fn, msg) {
  // Allow for deprecating things in the process of starting up.
  if (isUndefined(global.process)) {
    return function() {
      return exports.deprecate(fn, msg).apply(this, arguments);
    };
  }

  if (process.noDeprecation === true) {
    return fn;
  }

  var warned = false;
  function deprecated() {
    if (!warned) {
      if (process.throwDeprecation) {
        throw new Error(msg);
      } else if (process.traceDeprecation) {
        console.trace(msg);
      } else {
        console.error(msg);
      }
      warned = true;
    }
    return fn.apply(this, arguments);
  }

  return deprecated;
};


var debugs = {};
var debugEnviron;
exports.debuglog = function(set) {
  if (isUndefined(debugEnviron))
    debugEnviron = process.env.NODE_DEBUG || '';
  set = set.toUpperCase();
  if (!debugs[set]) {
    if (new RegExp('\\b' + set + '\\b', 'i').test(debugEnviron)) {
      var pid = process.pid;
      debugs[set] = function() {
        var msg = exports.format.apply(exports, arguments);
        console.error('%s %d: %s', set, pid, msg);
      };
    } else {
      debugs[set] = function() {};
    }
  }
  return debugs[set];
};


/**
 * Echos the value of a value. Trys to print the value out
 * in the best way possible given the different types.
 *
 * @param {Object} obj The object to print out.
 * @param {Object} opts Optional options object that alters the output.
 */
/* legacy: obj, showHidden, depth, colors*/
function inspect(obj, opts) {
  // default options
  var ctx = {
    seen: [],
    stylize: stylizeNoColor
  };
  // legacy...
  if (arguments.length >= 3) ctx.depth = arguments[2];
  if (arguments.length >= 4) ctx.colors = arguments[3];
  if (isBoolean(opts)) {
    // legacy...
    ctx.showHidden = opts;
  } else if (opts) {
    // got an "options" object
    exports._extend(ctx, opts);
  }
  // set default options
  if (isUndefined(ctx.showHidden)) ctx.showHidden = false;
  if (isUndefined(ctx.depth)) ctx.depth = 2;
  if (isUndefined(ctx.colors)) ctx.colors = false;
  if (isUndefined(ctx.customInspect)) ctx.customInspect = true;
  if (ctx.colors) ctx.stylize = stylizeWithColor;
  return formatValue(ctx, obj, ctx.depth);
}
exports.inspect = inspect;


// http://en.wikipedia.org/wiki/ANSI_escape_code#graphics
inspect.colors = {
  'bold' : [1, 22],
  'italic' : [3, 23],
  'underline' : [4, 24],
  'inverse' : [7, 27],
  'white' : [37, 39],
  'grey' : [90, 39],
  'black' : [30, 39],
  'blue' : [34, 39],
  'cyan' : [36, 39],
  'green' : [32, 39],
  'magenta' : [35, 39],
  'red' : [31, 39],
  'yellow' : [33, 39]
};

// Don't use 'blue' not visible on cmd.exe
inspect.styles = {
  'special': 'cyan',
  'number': 'yellow',
  'boolean': 'yellow',
  'undefined': 'grey',
  'null': 'bold',
  'string': 'green',
  'date': 'magenta',
  // "name": intentionally not styling
  'regexp': 'red'
};


function stylizeWithColor(str, styleType) {
  var style = inspect.styles[styleType];

  if (style) {
    return '\u001b[' + inspect.colors[style][0] + 'm' + str +
           '\u001b[' + inspect.colors[style][1] + 'm';
  } else {
    return str;
  }
}


function stylizeNoColor(str, styleType) {
  return str;
}


function arrayToHash(array) {
  var hash = {};

  array.forEach(function(val, idx) {
    hash[val] = true;
  });

  return hash;
}


function formatValue(ctx, value, recurseTimes) {
  // Provide a hook for user-specified inspect functions.
  // Check that value is an object with an inspect function on it
  if (ctx.customInspect &&
      value &&
      isFunction(value.inspect) &&
      // Filter out the util module, it's inspect function is special
      value.inspect !== exports.inspect &&
      // Also filter out any prototype objects using the circular check.
      !(value.constructor && value.constructor.prototype === value)) {
    var ret = value.inspect(recurseTimes, ctx);
    if (!isString(ret)) {
      ret = formatValue(ctx, ret, recurseTimes);
    }
    return ret;
  }

  // Primitive types cannot have properties
  var primitive = formatPrimitive(ctx, value);
  if (primitive) {
    return primitive;
  }

  // Look up the keys of the object.
  var keys = Object.keys(value);
  var visibleKeys = arrayToHash(keys);

  if (ctx.showHidden) {
    keys = Object.getOwnPropertyNames(value);
  }

  // IE doesn't make error fields non-enumerable
  // http://msdn.microsoft.com/en-us/library/ie/dww52sbt(v=vs.94).aspx
  if (isError(value)
      && (keys.indexOf('message') >= 0 || keys.indexOf('description') >= 0)) {
    return formatError(value);
  }

  // Some type of object without properties can be shortcutted.
  if (keys.length === 0) {
    if (isFunction(value)) {
      var name = value.name ? ': ' + value.name : '';
      return ctx.stylize('[Function' + name + ']', 'special');
    }
    if (isRegExp(value)) {
      return ctx.stylize(RegExp.prototype.toString.call(value), 'regexp');
    }
    if (isDate(value)) {
      return ctx.stylize(Date.prototype.toString.call(value), 'date');
    }
    if (isError(value)) {
      return formatError(value);
    }
  }

  var base = '', array = false, braces = ['{', '}'];

  // Make Array say that they are Array
  if (isArray(value)) {
    array = true;
    braces = ['[', ']'];
  }

  // Make functions say that they are functions
  if (isFunction(value)) {
    var n = value.name ? ': ' + value.name : '';
    base = ' [Function' + n + ']';
  }

  // Make RegExps say that they are RegExps
  if (isRegExp(value)) {
    base = ' ' + RegExp.prototype.toString.call(value);
  }

  // Make dates with properties first say the date
  if (isDate(value)) {
    base = ' ' + Date.prototype.toUTCString.call(value);
  }

  // Make error with message first say the error
  if (isError(value)) {
    base = ' ' + formatError(value);
  }

  if (keys.length === 0 && (!array || value.length == 0)) {
    return braces[0] + base + braces[1];
  }

  if (recurseTimes < 0) {
    if (isRegExp(value)) {
      return ctx.stylize(RegExp.prototype.toString.call(value), 'regexp');
    } else {
      return ctx.stylize('[Object]', 'special');
    }
  }

  ctx.seen.push(value);

  var output;
  if (array) {
    output = formatArray(ctx, value, recurseTimes, visibleKeys, keys);
  } else {
    output = keys.map(function(key) {
      return formatProperty(ctx, value, recurseTimes, visibleKeys, key, array);
    });
  }

  ctx.seen.pop();

  return reduceToSingleString(output, base, braces);
}


function formatPrimitive(ctx, value) {
  if (isUndefined(value))
    return ctx.stylize('undefined', 'undefined');
  if (isString(value)) {
    var simple = '\'' + JSON.stringify(value).replace(/^"|"$/g, '')
                                             .replace(/'/g, "\\'")
                                             .replace(/\\"/g, '"') + '\'';
    return ctx.stylize(simple, 'string');
  }
  if (isNumber(value))
    return ctx.stylize('' + value, 'number');
  if (isBoolean(value))
    return ctx.stylize('' + value, 'boolean');
  // For some reason typeof null is "object", so special case here.
  if (isNull(value))
    return ctx.stylize('null', 'null');
}


function formatError(value) {
  return '[' + Error.prototype.toString.call(value) + ']';
}


function formatArray(ctx, value, recurseTimes, visibleKeys, keys) {
  var output = [];
  for (var i = 0, l = value.length; i < l; ++i) {
    if (hasOwnProperty(value, String(i))) {
      output.push(formatProperty(ctx, value, recurseTimes, visibleKeys,
          String(i), true));
    } else {
      output.push('');
    }
  }
  keys.forEach(function(key) {
    if (!key.match(/^\d+$/)) {
      output.push(formatProperty(ctx, value, recurseTimes, visibleKeys,
          key, true));
    }
  });
  return output;
}


function formatProperty(ctx, value, recurseTimes, visibleKeys, key, array) {
  var name, str, desc;
  desc = Object.getOwnPropertyDescriptor(value, key) || { value: value[key] };
  if (desc.get) {
    if (desc.set) {
      str = ctx.stylize('[Getter/Setter]', 'special');
    } else {
      str = ctx.stylize('[Getter]', 'special');
    }
  } else {
    if (desc.set) {
      str = ctx.stylize('[Setter]', 'special');
    }
  }
  if (!hasOwnProperty(visibleKeys, key)) {
    name = '[' + key + ']';
  }
  if (!str) {
    if (ctx.seen.indexOf(desc.value) < 0) {
      if (isNull(recurseTimes)) {
        str = formatValue(ctx, desc.value, null);
      } else {
        str = formatValue(ctx, desc.value, recurseTimes - 1);
      }
      if (str.indexOf('\n') > -1) {
        if (array) {
          str = str.split('\n').map(function(line) {
            return '  ' + line;
          }).join('\n').substr(2);
        } else {
          str = '\n' + str.split('\n').map(function(line) {
            return '   ' + line;
          }).join('\n');
        }
      }
    } else {
      str = ctx.stylize('[Circular]', 'special');
    }
  }
  if (isUndefined(name)) {
    if (array && key.match(/^\d+$/)) {
      return str;
    }
    name = JSON.stringify('' + key);
    if (name.match(/^"([a-zA-Z_][a-zA-Z_0-9]*)"$/)) {
      name = name.substr(1, name.length - 2);
      name = ctx.stylize(name, 'name');
    } else {
      name = name.replace(/'/g, "\\'")
                 .replace(/\\"/g, '"')
                 .replace(/(^"|"$)/g, "'");
      name = ctx.stylize(name, 'string');
    }
  }

  return name + ': ' + str;
}


function reduceToSingleString(output, base, braces) {
  var numLinesEst = 0;
  var length = output.reduce(function(prev, cur) {
    numLinesEst++;
    if (cur.indexOf('\n') >= 0) numLinesEst++;
    return prev + cur.replace(/\u001b\[\d\d?m/g, '').length + 1;
  }, 0);

  if (length > 60) {
    return braces[0] +
           (base === '' ? '' : base + '\n ') +
           ' ' +
           output.join(',\n  ') +
           ' ' +
           braces[1];
  }

  return braces[0] + base + ' ' + output.join(', ') + ' ' + braces[1];
}


// NOTE: These type checking functions intentionally don't use `instanceof`
// because it is fragile and can be easily faked with `Object.create()`.
function isArray(ar) {
  return Array.isArray(ar);
}
exports.isArray = isArray;

function isBoolean(arg) {
  return typeof arg === 'boolean';
}
exports.isBoolean = isBoolean;

function isNull(arg) {
  return arg === null;
}
exports.isNull = isNull;

function isNullOrUndefined(arg) {
  return arg == null;
}
exports.isNullOrUndefined = isNullOrUndefined;

function isNumber(arg) {
  return typeof arg === 'number';
}
exports.isNumber = isNumber;

function isString(arg) {
  return typeof arg === 'string';
}
exports.isString = isString;

function isSymbol(arg) {
  return typeof arg === 'symbol';
}
exports.isSymbol = isSymbol;

function isUndefined(arg) {
  return arg === void 0;
}
exports.isUndefined = isUndefined;

function isRegExp(re) {
  return isObject(re) && objectToString(re) === '[object RegExp]';
}
exports.isRegExp = isRegExp;

function isObject(arg) {
  return typeof arg === 'object' && arg !== null;
}
exports.isObject = isObject;

function isDate(d) {
  return isObject(d) && objectToString(d) === '[object Date]';
}
exports.isDate = isDate;

function isError(e) {
  return isObject(e) &&
      (objectToString(e) === '[object Error]' || e instanceof Error);
}
exports.isError = isError;

function isFunction(arg) {
  return typeof arg === 'function';
}
exports.isFunction = isFunction;

function isPrimitive(arg) {
  return arg === null ||
         typeof arg === 'boolean' ||
         typeof arg === 'number' ||
         typeof arg === 'string' ||
         typeof arg === 'symbol' ||  // ES6 symbol
         typeof arg === 'undefined';
}
exports.isPrimitive = isPrimitive;

exports.isBuffer = require('./support/isBuffer');

function objectToString(o) {
  return Object.prototype.toString.call(o);
}


function pad(n) {
  return n < 10 ? '0' + n.toString(10) : n.toString(10);
}


var months = ['Jan', 'Feb', 'Mar', 'Apr', 'May', 'Jun', 'Jul', 'Aug', 'Sep',
              'Oct', 'Nov', 'Dec'];

// 26 Feb 16:19:34
function timestamp() {
  var d = new Date();
  var time = [pad(d.getHours()),
              pad(d.getMinutes()),
              pad(d.getSeconds())].join(':');
  return [d.getDate(), months[d.getMonth()], time].join(' ');
}


// log is just a thin wrapper to console.log that prepends a timestamp
exports.log = function() {
  console.log('%s - %s', timestamp(), exports.format.apply(exports, arguments));
};


/**
 * Inherit the prototype methods from one constructor into another.
 *
 * The Function.prototype.inherits from lang.js rewritten as a standalone
 * function (not on Function.prototype). NOTE: If this file is to be loaded
 * during bootstrapping this function needs to be rewritten using some native
 * functions as prototype setup using normal JavaScript does not work as
 * expected during bootstrapping (see mirror.js in r114903).
 *
 * @param {function} ctor Constructor function which needs to inherit the
 *     prototype.
 * @param {function} superCtor Constructor function to inherit prototype from.
 */
exports.inherits = require('inherits');

exports._extend = function(origin, add) {
  // Don't do anything if add isn't an object
  if (!add || !isObject(add)) return origin;

  var keys = Object.keys(add);
  var i = keys.length;
  while (i--) {
    origin[keys[i]] = add[keys[i]];
  }
  return origin;
};

function hasOwnProperty(obj, prop) {
  return Object.prototype.hasOwnProperty.call(obj, prop);
}

}).call(this,require('_process'),typeof global !== "undefined" ? global : typeof self !== "undefined" ? self : typeof window !== "undefined" ? window : {})
},{"./support/isBuffer":"/usr/local/lib/node_modules/watchify/node_modules/browserify/node_modules/util/support/isBufferBrowser.js","_process":"/usr/local/lib/node_modules/watchify/node_modules/browserify/node_modules/process/browser.js","inherits":"/usr/local/lib/node_modules/watchify/node_modules/browserify/node_modules/inherits/inherits_browser.js"}]},{},["/Users/Romain/Documents/Programmation/Projets_Web/AntColony/src/index.js"])
//# sourceMappingURL=data:application/json;base64,eyJ2ZXJzaW9uIjozLCJzb3VyY2VzIjpbIi91c3IvbG9jYWwvbGliL25vZGVfbW9kdWxlcy93YXRjaGlmeS9ub2RlX21vZHVsZXMvYnJvd3NlcmlmeS9ub2RlX21vZHVsZXMvYnJvd3Nlci1wYWNrL19wcmVsdWRlLmpzIiwiL1VzZXJzL1JvbWFpbi9Eb2N1bWVudHMvUHJvZ3JhbW1hdGlvbi9Qcm9qZXRzX1dlYi9BbnRDb2xvbnkvbm9kZV9tb2R1bGVzL2RlbGF1bmF5LXRyaWFuZ3VsYXRlL25vZGVfbW9kdWxlcy9pbmNyZW1lbnRhbC1jb252ZXgtaHVsbC9pY2guanMiLCIvVXNlcnMvUm9tYWluL0RvY3VtZW50cy9Qcm9ncmFtbWF0aW9uL1Byb2pldHNfV2ViL0FudENvbG9ueS9ub2RlX21vZHVsZXMvZGVsYXVuYXktdHJpYW5ndWxhdGUvbm9kZV9tb2R1bGVzL2luY3JlbWVudGFsLWNvbnZleC1odWxsL25vZGVfbW9kdWxlcy9yb2J1c3Qtb3JpZW50YXRpb24vbm9kZV9tb2R1bGVzL3JvYnVzdC1zY2FsZS9ub2RlX21vZHVsZXMvdHdvLXN1bS90d28tc3VtLmpzIiwiL1VzZXJzL1JvbWFpbi9Eb2N1bWVudHMvUHJvZ3JhbW1hdGlvbi9Qcm9qZXRzX1dlYi9BbnRDb2xvbnkvbm9kZV9tb2R1bGVzL2RlbGF1bmF5LXRyaWFuZ3VsYXRlL25vZGVfbW9kdWxlcy9pbmNyZW1lbnRhbC1jb252ZXgtaHVsbC9ub2RlX21vZHVsZXMvcm9idXN0LW9yaWVudGF0aW9uL25vZGVfbW9kdWxlcy9yb2J1c3Qtc2NhbGUvcm9idXN0LXNjYWxlLmpzIiwiL1VzZXJzL1JvbWFpbi9Eb2N1bWVudHMvUHJvZ3JhbW1hdGlvbi9Qcm9qZXRzX1dlYi9BbnRDb2xvbnkvbm9kZV9tb2R1bGVzL2RlbGF1bmF5LXRyaWFuZ3VsYXRlL25vZGVfbW9kdWxlcy9pbmNyZW1lbnRhbC1jb252ZXgtaHVsbC9ub2RlX21vZHVsZXMvcm9idXN0LW9yaWVudGF0aW9uL25vZGVfbW9kdWxlcy9yb2J1c3Qtc3VidHJhY3Qvcm9idXN0LWRpZmYuanMiLCIvVXNlcnMvUm9tYWluL0RvY3VtZW50cy9Qcm9ncmFtbWF0aW9uL1Byb2pldHNfV2ViL0FudENvbG9ueS9ub2RlX21vZHVsZXMvZGVsYXVuYXktdHJpYW5ndWxhdGUvbm9kZV9tb2R1bGVzL2luY3JlbWVudGFsLWNvbnZleC1odWxsL25vZGVfbW9kdWxlcy9yb2J1c3Qtb3JpZW50YXRpb24vbm9kZV9tb2R1bGVzL3JvYnVzdC1zdW0vcm9idXN0LXN1bS5qcyIsIi9Vc2Vycy9Sb21haW4vRG9jdW1lbnRzL1Byb2dyYW1tYXRpb24vUHJvamV0c19XZWIvQW50Q29sb255L25vZGVfbW9kdWxlcy9kZWxhdW5heS10cmlhbmd1bGF0ZS9ub2RlX21vZHVsZXMvaW5jcmVtZW50YWwtY29udmV4LWh1bGwvbm9kZV9tb2R1bGVzL3JvYnVzdC1vcmllbnRhdGlvbi9ub2RlX21vZHVsZXMvdHdvLXByb2R1Y3QvdHdvLXByb2R1Y3QuanMiLCIvVXNlcnMvUm9tYWluL0RvY3VtZW50cy9Qcm9ncmFtbWF0aW9uL1Byb2pldHNfV2ViL0FudENvbG9ueS9ub2RlX21vZHVsZXMvZGVsYXVuYXktdHJpYW5ndWxhdGUvbm9kZV9tb2R1bGVzL2luY3JlbWVudGFsLWNvbnZleC1odWxsL25vZGVfbW9kdWxlcy9yb2J1c3Qtb3JpZW50YXRpb24vb3JpZW50YXRpb24uanMiLCIvVXNlcnMvUm9tYWluL0RvY3VtZW50cy9Qcm9ncmFtbWF0aW9uL1Byb2pldHNfV2ViL0FudENvbG9ueS9ub2RlX21vZHVsZXMvZGVsYXVuYXktdHJpYW5ndWxhdGUvbm9kZV9tb2R1bGVzL2luY3JlbWVudGFsLWNvbnZleC1odWxsL25vZGVfbW9kdWxlcy9zaW1wbGljaWFsLWNvbXBsZXgvbm9kZV9tb2R1bGVzL2JpdC10d2lkZGxlL3R3aWRkbGUuanMiLCIvVXNlcnMvUm9tYWluL0RvY3VtZW50cy9Qcm9ncmFtbWF0aW9uL1Byb2pldHNfV2ViL0FudENvbG9ueS9ub2RlX21vZHVsZXMvZGVsYXVuYXktdHJpYW5ndWxhdGUvbm9kZV9tb2R1bGVzL2luY3JlbWVudGFsLWNvbnZleC1odWxsL25vZGVfbW9kdWxlcy9zaW1wbGljaWFsLWNvbXBsZXgvbm9kZV9tb2R1bGVzL3VuaW9uLWZpbmQvaW5kZXguanMiLCIvVXNlcnMvUm9tYWluL0RvY3VtZW50cy9Qcm9ncmFtbWF0aW9uL1Byb2pldHNfV2ViL0FudENvbG9ueS9ub2RlX21vZHVsZXMvZGVsYXVuYXktdHJpYW5ndWxhdGUvbm9kZV9tb2R1bGVzL2luY3JlbWVudGFsLWNvbnZleC1odWxsL25vZGVfbW9kdWxlcy9zaW1wbGljaWFsLWNvbXBsZXgvdG9wb2xvZ3kuanMiLCIvVXNlcnMvUm9tYWluL0RvY3VtZW50cy9Qcm9ncmFtbWF0aW9uL1Byb2pldHNfV2ViL0FudENvbG9ueS9ub2RlX21vZHVsZXMvZGVsYXVuYXktdHJpYW5ndWxhdGUvbm9kZV9tb2R1bGVzL3VuaXEvdW5pcS5qcyIsIi9Vc2Vycy9Sb21haW4vRG9jdW1lbnRzL1Byb2dyYW1tYXRpb24vUHJvamV0c19XZWIvQW50Q29sb255L25vZGVfbW9kdWxlcy9kZWxhdW5heS10cmlhbmd1bGF0ZS90cmlhbmd1bGF0ZS5qcyIsIi9Vc2Vycy9Sb21haW4vRG9jdW1lbnRzL1Byb2dyYW1tYXRpb24vUHJvamV0c19XZWIvQW50Q29sb255L25vZGVfbW9kdWxlcy9nYW1lLXNoZWxsL2xpYi9ocnRpbWUtcG9seWZpbGwuanMiLCIvVXNlcnMvUm9tYWluL0RvY3VtZW50cy9Qcm9ncmFtbWF0aW9uL1Byb2pldHNfV2ViL0FudENvbG9ueS9ub2RlX21vZHVsZXMvZ2FtZS1zaGVsbC9saWIvbW91c2V3aGVlbC1wb2x5ZmlsbC5qcyIsIi9Vc2Vycy9Sb21haW4vRG9jdW1lbnRzL1Byb2dyYW1tYXRpb24vUHJvamV0c19XZWIvQW50Q29sb255L25vZGVfbW9kdWxlcy9nYW1lLXNoZWxsL2xpYi9yYWYtcG9seWZpbGwuanMiLCIvVXNlcnMvUm9tYWluL0RvY3VtZW50cy9Qcm9ncmFtbWF0aW9uL1Byb2pldHNfV2ViL0FudENvbG9ueS9ub2RlX21vZHVsZXMvZ2FtZS1zaGVsbC9ub2RlX21vZHVsZXMvYmluYXJ5LXNlYXJjaC1ib3VuZHMvc2VhcmNoLWJvdW5kcy5qcyIsIi9Vc2Vycy9Sb21haW4vRG9jdW1lbnRzL1Byb2dyYW1tYXRpb24vUHJvamV0c19XZWIvQW50Q29sb255L25vZGVfbW9kdWxlcy9nYW1lLXNoZWxsL25vZGVfbW9kdWxlcy9kb21yZWFkeS9yZWFkeS5qcyIsIi9Vc2Vycy9Sb21haW4vRG9jdW1lbnRzL1Byb2dyYW1tYXRpb24vUHJvamV0c19XZWIvQW50Q29sb255L25vZGVfbW9kdWxlcy9nYW1lLXNoZWxsL25vZGVfbW9kdWxlcy9pbnZlcnQtaGFzaC9pbnZlcnQuanMiLCIvVXNlcnMvUm9tYWluL0RvY3VtZW50cy9Qcm9ncmFtbWF0aW9uL1Byb2pldHNfV2ViL0FudENvbG9ueS9ub2RlX21vZHVsZXMvZ2FtZS1zaGVsbC9ub2RlX21vZHVsZXMvaW90YS1hcnJheS9pb3RhLmpzIiwiL1VzZXJzL1JvbWFpbi9Eb2N1bWVudHMvUHJvZ3JhbW1hdGlvbi9Qcm9qZXRzX1dlYi9BbnRDb2xvbnkvbm9kZV9tb2R1bGVzL2dhbWUtc2hlbGwvbm9kZV9tb2R1bGVzL3VuaXEvdW5pcS5qcyIsIi9Vc2Vycy9Sb21haW4vRG9jdW1lbnRzL1Byb2dyYW1tYXRpb24vUHJvamV0c19XZWIvQW50Q29sb255L25vZGVfbW9kdWxlcy9nYW1lLXNoZWxsL25vZGVfbW9kdWxlcy92a2V5L2luZGV4LmpzIiwiL1VzZXJzL1JvbWFpbi9Eb2N1bWVudHMvUHJvZ3JhbW1hdGlvbi9Qcm9qZXRzX1dlYi9BbnRDb2xvbnkvbm9kZV9tb2R1bGVzL2dhbWUtc2hlbGwvc2hlbGwuanMiLCIvVXNlcnMvUm9tYWluL0RvY3VtZW50cy9Qcm9ncmFtbWF0aW9uL1Byb2pldHNfV2ViL0FudENvbG9ueS9ub2RlX21vZHVsZXMvaGFybW9ueS1jb2xsZWN0aW9ucy9oYXJtb255LWNvbGxlY3Rpb25zLmpzIiwiL1VzZXJzL1JvbWFpbi9Eb2N1bWVudHMvUHJvZ3JhbW1hdGlvbi9Qcm9qZXRzX1dlYi9BbnRDb2xvbnkvbm9kZV9tb2R1bGVzL3BhcnNlLXN2Zy1wYXRoL2luZGV4LmpzIiwiL1VzZXJzL1JvbWFpbi9Eb2N1bWVudHMvUHJvZ3JhbW1hdGlvbi9Qcm9qZXRzX1dlYi9BbnRDb2xvbnkvc3JjL2luZGV4LmpzIiwiL3Vzci9sb2NhbC9saWIvbm9kZV9tb2R1bGVzL3dhdGNoaWZ5L25vZGVfbW9kdWxlcy9icm93c2VyaWZ5L25vZGVfbW9kdWxlcy9ldmVudHMvZXZlbnRzLmpzIiwiL3Vzci9sb2NhbC9saWIvbm9kZV9tb2R1bGVzL3dhdGNoaWZ5L25vZGVfbW9kdWxlcy9icm93c2VyaWZ5L25vZGVfbW9kdWxlcy9pbmhlcml0cy9pbmhlcml0c19icm93c2VyLmpzIiwiL3Vzci9sb2NhbC9saWIvbm9kZV9tb2R1bGVzL3dhdGNoaWZ5L25vZGVfbW9kdWxlcy9icm93c2VyaWZ5L25vZGVfbW9kdWxlcy9wcm9jZXNzL2Jyb3dzZXIuanMiLCIvdXNyL2xvY2FsL2xpYi9ub2RlX21vZHVsZXMvd2F0Y2hpZnkvbm9kZV9tb2R1bGVzL2Jyb3dzZXJpZnkvbm9kZV9tb2R1bGVzL3V0aWwvc3VwcG9ydC9pc0J1ZmZlckJyb3dzZXIuanMiLCIvdXNyL2xvY2FsL2xpYi9ub2RlX21vZHVsZXMvd2F0Y2hpZnkvbm9kZV9tb2R1bGVzL2Jyb3dzZXJpZnkvbm9kZV9tb2R1bGVzL3V0aWwvdXRpbC5qcyJdLCJuYW1lcyI6W10sIm1hcHBpbmdzIjoiQUFBQTtBQ0FBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7O0FDN2JBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7O0FDaEJBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7O0FDakRBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTs7QUMzSkE7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBOztBQzNKQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7O0FDaENBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBOztBQzdMQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTs7QUM1TUE7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7O0FDekRBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBOztBQ3RWQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTs7QUN6REE7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBOztBQzlKQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7O0FDWEE7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTs7QUMxREE7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTs7QUM1QkE7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7O0FDNURBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7O0FDNUJBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBOztBQ1pBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7O0FDVkE7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7O0FDekRBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7O0FDeElBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7O0FDOXNCQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTs7QUMzeEJBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7O0FDdkRBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTs7QUNqVEE7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7O0FDL1NBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTs7QUN2QkE7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7O0FDL0RBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTs7QUNMQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQSIsImZpbGUiOiJnZW5lcmF0ZWQuanMiLCJzb3VyY2VSb290IjoiIiwic291cmNlc0NvbnRlbnQiOlsiKGZ1bmN0aW9uIGUodCxuLHIpe2Z1bmN0aW9uIHMobyx1KXtpZighbltvXSl7aWYoIXRbb10pe3ZhciBhPXR5cGVvZiByZXF1aXJlPT1cImZ1bmN0aW9uXCImJnJlcXVpcmU7aWYoIXUmJmEpcmV0dXJuIGEobywhMCk7aWYoaSlyZXR1cm4gaShvLCEwKTt2YXIgZj1uZXcgRXJyb3IoXCJDYW5ub3QgZmluZCBtb2R1bGUgJ1wiK28rXCInXCIpO3Rocm93IGYuY29kZT1cIk1PRFVMRV9OT1RfRk9VTkRcIixmfXZhciBsPW5bb109e2V4cG9ydHM6e319O3Rbb11bMF0uY2FsbChsLmV4cG9ydHMsZnVuY3Rpb24oZSl7dmFyIG49dFtvXVsxXVtlXTtyZXR1cm4gcyhuP246ZSl9LGwsbC5leHBvcnRzLGUsdCxuLHIpfXJldHVybiBuW29dLmV4cG9ydHN9dmFyIGk9dHlwZW9mIHJlcXVpcmU9PVwiZnVuY3Rpb25cIiYmcmVxdWlyZTtmb3IodmFyIG89MDtvPHIubGVuZ3RoO28rKylzKHJbb10pO3JldHVybiBzfSkiLCJcInVzZSBzdHJpY3RcIlxuXG4vL0hpZ2ggbGV2ZWwgaWRlYTpcbi8vIDEuIFVzZSBDbGFya3NvbidzIGluY3JlbWVudGFsIGNvbnN0cnVjdGlvbiB0byBmaW5kIGNvbnZleCBodWxsXG4vLyAyLiBQb2ludCBsb2NhdGlvbiBpbiB0cmlhbmd1bGF0aW9uIGJ5IGp1bXAgYW5kIHdhbGtcblxubW9kdWxlLmV4cG9ydHMgPSBpbmNyZW1lbnRhbENvbnZleEh1bGxcblxudmFyIG9yaWVudCA9IHJlcXVpcmUoXCJyb2J1c3Qtb3JpZW50YXRpb25cIilcbnZhciBjb21wYXJlQ2VsbCA9IHJlcXVpcmUoXCJzaW1wbGljaWFsLWNvbXBsZXhcIikuY29tcGFyZUNlbGxzXG5cbmZ1bmN0aW9uIGNvbXBhcmVJbnQoYSwgYikge1xuICByZXR1cm4gYSAtIGJcbn1cblxuZnVuY3Rpb24gU2ltcGxleCh2ZXJ0aWNlcywgYWRqYWNlbnQsIGJvdW5kYXJ5KSB7XG4gIHRoaXMudmVydGljZXMgPSB2ZXJ0aWNlc1xuICB0aGlzLmFkamFjZW50ID0gYWRqYWNlbnRcbiAgdGhpcy5ib3VuZGFyeSA9IGJvdW5kYXJ5XG4gIHRoaXMubGFzdFZpc2l0ZWQgPSAtMVxufVxuXG5TaW1wbGV4LnByb3RvdHlwZS5mbGlwID0gZnVuY3Rpb24oKSB7XG4gIHZhciB0ID0gdGhpcy52ZXJ0aWNlc1swXVxuICB0aGlzLnZlcnRpY2VzWzBdID0gdGhpcy52ZXJ0aWNlc1sxXVxuICB0aGlzLnZlcnRpY2VzWzFdID0gdFxuICB2YXIgdSA9IHRoaXMuYWRqYWNlbnRbMF1cbiAgdGhpcy5hZGphY2VudFswXSA9IHRoaXMuYWRqYWNlbnRbMV1cbiAgdGhpcy5hZGphY2VudFsxXSA9IHVcbn1cblxuZnVuY3Rpb24gR2x1ZUZhY2V0KHZlcnRpY2VzLCBjZWxsLCBpbmRleCkge1xuICB0aGlzLnZlcnRpY2VzID0gdmVydGljZXNcbiAgdGhpcy5jZWxsID0gY2VsbFxuICB0aGlzLmluZGV4ID0gaW5kZXhcbn1cblxuZnVuY3Rpb24gY29tcGFyZUdsdWUoYSwgYikge1xuICByZXR1cm4gY29tcGFyZUNlbGwoYS52ZXJ0aWNlcywgYi52ZXJ0aWNlcylcbn1cblxuZnVuY3Rpb24gYmFrZU9yaWVudChkKSB7XG4gIHZhciBjb2RlID0gW1wiZnVuY3Rpb24gb3JpZW50KCl7dmFyIHR1cGxlPXRoaXMudHVwbGU7cmV0dXJuIHRlc3QoXCJdXG4gIGZvcih2YXIgaT0wOyBpPD1kOyArK2kpIHtcbiAgICBpZihpID4gMCkge1xuICAgICAgY29kZS5wdXNoKFwiLFwiKVxuICAgIH1cbiAgICBjb2RlLnB1c2goXCJ0dXBsZVtcIiwgaSwgXCJdXCIpXG4gIH1cbiAgY29kZS5wdXNoKFwiKX1yZXR1cm4gb3JpZW50XCIpXG4gIHZhciBwcm9jID0gbmV3IEZ1bmN0aW9uKFwidGVzdFwiLCBjb2RlLmpvaW4oXCJcIikpXG4gIHZhciB0ZXN0ID0gb3JpZW50W2QrMV1cbiAgaWYoIXRlc3QpIHtcbiAgICB0ZXN0ID0gb3JpZW50XG4gIH1cbiAgcmV0dXJuIHByb2ModGVzdClcbn1cblxudmFyIEJBS0VEID0gW11cblxuZnVuY3Rpb24gVHJpYW5ndWxhdGlvbihkaW1lbnNpb24sIHZlcnRpY2VzLCBzaW1wbGljZXMpIHtcbiAgdGhpcy5kaW1lbnNpb24gPSBkaW1lbnNpb25cbiAgdGhpcy52ZXJ0aWNlcyA9IHZlcnRpY2VzXG4gIHRoaXMuc2ltcGxpY2VzID0gc2ltcGxpY2VzXG4gIHRoaXMuaW50ZXJpb3IgPSBzaW1wbGljZXMuZmlsdGVyKGZ1bmN0aW9uKGMpIHtcbiAgICByZXR1cm4gIWMuYm91bmRhcnlcbiAgfSlcblxuICB0aGlzLnR1cGxlID0gbmV3IEFycmF5KGRpbWVuc2lvbisxKVxuICBmb3IodmFyIGk9MDsgaTw9ZGltZW5zaW9uOyArK2kpIHtcbiAgICB0aGlzLnR1cGxlW2ldID0gdGhpcy52ZXJ0aWNlc1tpXVxuICB9XG5cbiAgdmFyIG8gPSBCQUtFRFtkaW1lbnNpb25dXG4gIGlmKCFvKSB7XG4gICAgbyA9IEJBS0VEW2RpbWVuc2lvbl0gPSBiYWtlT3JpZW50KGRpbWVuc2lvbilcbiAgfVxuICB0aGlzLm9yaWVudCA9IG9cbn1cblxudmFyIHByb3RvID0gVHJpYW5ndWxhdGlvbi5wcm90b3R5cGVcblxuLy9EZWdlbmVyYXRlIHNpdHVhdGlvbiB3aGVyZSB3ZSBhcmUgb24gYm91bmRhcnksIGJ1dCBjb3BsYW5hciB0byBmYWNlXG5wcm90by5oYW5kbGVCb3VuZGFyeURlZ2VuZXJhY3kgPSBmdW5jdGlvbihjZWxsLCBwb2ludCkge1xuICB2YXIgZCA9IHRoaXMuZGltZW5zaW9uXG4gIHZhciBuID0gdGhpcy52ZXJ0aWNlcy5sZW5ndGggLSAxXG4gIHZhciB0dXBsZSA9IHRoaXMudHVwbGVcbiAgdmFyIHZlcnRzID0gdGhpcy52ZXJ0aWNlc1xuXG4gIC8vRHVtYiBzb2x1dGlvbjogSnVzdCBkbyBkZnMgZnJvbSBib3VuZGFyeSBjZWxsIHVudGlsIHdlIGZpbmQgYW55IHBlYWssIG9yIHRlcm1pbmF0ZVxuICB2YXIgdG9WaXNpdCA9IFsgY2VsbCBdXG4gIGNlbGwubGFzdFZpc2l0ZWQgPSAtblxuICB3aGlsZSh0b1Zpc2l0Lmxlbmd0aCA+IDApIHtcbiAgICBjZWxsID0gdG9WaXNpdC5wb3AoKVxuICAgIHZhciBjZWxsVmVydHMgPSBjZWxsLnZlcnRpY2VzXG4gICAgdmFyIGNlbGxBZGogPSBjZWxsLmFkamFjZW50XG4gICAgZm9yKHZhciBpPTA7IGk8PWQ7ICsraSkge1xuICAgICAgdmFyIG5laWdoYm9yID0gY2VsbEFkaltpXVxuICAgICAgaWYoIW5laWdoYm9yLmJvdW5kYXJ5IHx8IG5laWdoYm9yLmxhc3RWaXNpdGVkIDw9IC1uKSB7XG4gICAgICAgIGNvbnRpbnVlXG4gICAgICB9XG4gICAgICB2YXIgbnYgPSBuZWlnaGJvci52ZXJ0aWNlc1xuICAgICAgZm9yKHZhciBqPTA7IGo8PWQ7ICsraikge1xuICAgICAgICB2YXIgdnYgPSBudltqXVxuICAgICAgICBpZih2diA8IDApIHtcbiAgICAgICAgICB0dXBsZVtqXSA9IHBvaW50XG4gICAgICAgIH0gZWxzZSB7XG4gICAgICAgICAgdHVwbGVbal0gPSB2ZXJ0c1t2dl1cbiAgICAgICAgfVxuICAgICAgfVxuICAgICAgdmFyIG8gPSB0aGlzLm9yaWVudCgpXG4gICAgICBpZihvID4gMCkge1xuICAgICAgICByZXR1cm4gbmVpZ2hib3JcbiAgICAgIH1cbiAgICAgIG5laWdoYm9yLmxhc3RWaXNpdGVkID0gLW5cbiAgICAgIGlmKG8gPT09IDApIHtcbiAgICAgICAgdG9WaXNpdC5wdXNoKG5laWdoYm9yKVxuICAgICAgfVxuICAgIH1cbiAgfVxuICByZXR1cm4gbnVsbFxufVxuXG5wcm90by53YWxrID0gZnVuY3Rpb24ocG9pbnQsIHJhbmRvbSkge1xuICAvL0FsaWFzIGxvY2FsIHByb3BlcnRpZXNcbiAgdmFyIG4gPSB0aGlzLnZlcnRpY2VzLmxlbmd0aCAtIDFcbiAgdmFyIGQgPSB0aGlzLmRpbWVuc2lvblxuICB2YXIgdmVydHMgPSB0aGlzLnZlcnRpY2VzXG4gIHZhciB0dXBsZSA9IHRoaXMudHVwbGVcblxuICAvL0NvbXB1dGUgaW5pdGlhbCBqdW1wIGNlbGxcbiAgdmFyIGluaXRJbmRleCA9IHJhbmRvbSA/ICh0aGlzLmludGVyaW9yLmxlbmd0aCAqIE1hdGgucmFuZG9tKCkpfDAgOiAodGhpcy5pbnRlcmlvci5sZW5ndGgtMSlcbiAgdmFyIGNlbGwgPSB0aGlzLmludGVyaW9yWyBpbml0SW5kZXggXVxuXG4gIC8vU3RhcnQgd2Fsa2luZ1xub3V0ZXJMb29wOlxuICB3aGlsZSghY2VsbC5ib3VuZGFyeSkge1xuICAgIHZhciBjZWxsVmVydHMgPSBjZWxsLnZlcnRpY2VzXG4gICAgdmFyIGNlbGxBZGogPSBjZWxsLmFkamFjZW50XG5cbiAgICBmb3IodmFyIGk9MDsgaTw9ZDsgKytpKSB7XG4gICAgICB0dXBsZVtpXSA9IHZlcnRzW2NlbGxWZXJ0c1tpXV1cbiAgICB9XG4gICAgY2VsbC5sYXN0VmlzaXRlZCA9IG5cblxuICAgIC8vRmluZCBmYXJ0aGVzdCBhZGphY2VudCBjZWxsXG4gICAgZm9yKHZhciBpPTA7IGk8PWQ7ICsraSkge1xuICAgICAgdmFyIG5laWdoYm9yID0gY2VsbEFkaltpXVxuICAgICAgaWYobmVpZ2hib3IubGFzdFZpc2l0ZWQgPj0gbikge1xuICAgICAgICBjb250aW51ZVxuICAgICAgfVxuICAgICAgdmFyIHByZXYgPSB0dXBsZVtpXVxuICAgICAgdHVwbGVbaV0gPSBwb2ludFxuICAgICAgdmFyIG8gPSB0aGlzLm9yaWVudCgpXG4gICAgICB0dXBsZVtpXSA9IHByZXZcbiAgICAgIGlmKG8gPCAwKSB7XG4gICAgICAgIGNlbGwgPSBuZWlnaGJvclxuICAgICAgICBjb250aW51ZSBvdXRlckxvb3BcbiAgICAgIH0gZWxzZSB7XG4gICAgICAgIGlmKCFuZWlnaGJvci5ib3VuZGFyeSkge1xuICAgICAgICAgIG5laWdoYm9yLmxhc3RWaXNpdGVkID0gblxuICAgICAgICB9IGVsc2Uge1xuICAgICAgICAgIG5laWdoYm9yLmxhc3RWaXNpdGVkID0gLW5cbiAgICAgICAgfVxuICAgICAgfVxuICAgIH1cbiAgICByZXR1cm5cbiAgfVxuXG4gIHJldHVybiBjZWxsXG59XG5cbnByb3RvLmFkZFBlYWtzID0gZnVuY3Rpb24ocG9pbnQsIGNlbGwpIHtcbiAgdmFyIG4gPSB0aGlzLnZlcnRpY2VzLmxlbmd0aCAtIDFcbiAgdmFyIGQgPSB0aGlzLmRpbWVuc2lvblxuICB2YXIgdmVydHMgPSB0aGlzLnZlcnRpY2VzXG4gIHZhciB0dXBsZSA9IHRoaXMudHVwbGVcbiAgdmFyIGludGVyaW9yID0gdGhpcy5pbnRlcmlvclxuICB2YXIgc2ltcGxpY2VzID0gdGhpcy5zaW1wbGljZXNcblxuICAvL1dhbGtpbmcgZmluaXNoZWQgYXQgYm91bmRhcnksIHRpbWUgdG8gYWRkIHBlYWtzXG4gIHZhciB0b3Zpc2l0ID0gWyBjZWxsIF1cblxuICAvL1N0cmV0Y2ggaW5pdGlhbCBib3VuZGFyeSBjZWxsIGludG8gYSBwZWFrXG4gIGNlbGwubGFzdFZpc2l0ZWQgPSBuXG4gIGNlbGwudmVydGljZXNbY2VsbC52ZXJ0aWNlcy5pbmRleE9mKC0xKV0gPSBuXG4gIGNlbGwuYm91bmRhcnkgPSBmYWxzZVxuICBpbnRlcmlvci5wdXNoKGNlbGwpXG5cbiAgLy9SZWNvcmQgYSBsaXN0IG9mIGFsbCBuZXcgYm91bmRhcmllcyBjcmVhdGVkIGJ5IGFkZGVkIHBlYWtzIHNvIHdlIGNhbiBnbHVlIHRoZW0gdG9nZXRoZXIgd2hlbiB3ZSBhcmUgYWxsIGRvbmVcbiAgdmFyIGdsdWVGYWNldHMgPSBbXVxuXG4gIC8vRG8gYSB0cmF2ZXJzYWwgb2YgdGhlIGJvdW5kYXJ5IHdhbGtpbmcgb3V0d2FyZCBmcm9tIHN0YXJ0aW5nIHBlYWtcbiAgd2hpbGUodG92aXNpdC5sZW5ndGggPiAwKSB7XG4gICAgLy9Qb3Agb2ZmIHBlYWsgYW5kIHdhbGsgb3ZlciBhZGphY2VudCBjZWxsc1xuICAgIHZhciBjZWxsID0gdG92aXNpdC5wb3AoKVxuICAgIHZhciBjZWxsVmVydHMgPSBjZWxsLnZlcnRpY2VzXG4gICAgdmFyIGNlbGxBZGogPSBjZWxsLmFkamFjZW50XG4gICAgdmFyIGluZGV4T2ZOID0gY2VsbFZlcnRzLmluZGV4T2YobilcbiAgICBpZihpbmRleE9mTiA8IDApIHtcbiAgICAgIGNvbnRpbnVlXG4gICAgfVxuXG4gICAgZm9yKHZhciBpPTA7IGk8PWQ7ICsraSkge1xuICAgICAgaWYoaSA9PT0gaW5kZXhPZk4pIHtcbiAgICAgICAgY29udGludWVcbiAgICAgIH1cblxuICAgICAgLy9Gb3IgZWFjaCBib3VuZGFyeSBuZWlnaGJvciBvZiB0aGUgY2VsbFxuICAgICAgdmFyIG5laWdoYm9yID0gY2VsbEFkaltpXVxuICAgICAgaWYoIW5laWdoYm9yLmJvdW5kYXJ5IHx8IG5laWdoYm9yLmxhc3RWaXNpdGVkID49IG4pIHtcbiAgICAgICAgY29udGludWVcbiAgICAgIH1cblxuICAgICAgdmFyIG52ID0gbmVpZ2hib3IudmVydGljZXNcblxuICAgICAgLy9UZXN0IGlmIG5laWdoYm9yIGlzIGEgcGVha1xuICAgICAgaWYobmVpZ2hib3IubGFzdFZpc2l0ZWQgIT09IC1uKSB7ICAgICAgXG4gICAgICAgIC8vQ29tcHV0ZSBvcmllbnRhdGlvbiBvZiBwIHJlbGF0aXZlIHRvIGVhY2ggYm91bmRhcnkgcGVha1xuICAgICAgICB2YXIgaW5kZXhPZk5lZzEgPSAwXG4gICAgICAgIGZvcih2YXIgaj0wOyBqPD1kOyArK2opIHtcbiAgICAgICAgICBpZihudltqXSA8IDApIHtcbiAgICAgICAgICAgIGluZGV4T2ZOZWcxID0galxuICAgICAgICAgICAgdHVwbGVbal0gPSBwb2ludFxuICAgICAgICAgIH0gZWxzZSB7XG4gICAgICAgICAgICB0dXBsZVtqXSA9IHZlcnRzW252W2pdXVxuICAgICAgICAgIH1cbiAgICAgICAgfVxuICAgICAgICB2YXIgbyA9IHRoaXMub3JpZW50KClcblxuICAgICAgICAvL1Rlc3QgaWYgbmVpZ2hib3IgY2VsbCBpcyBhbHNvIGEgcGVha1xuICAgICAgICBpZihvID4gMCkge1xuICAgICAgICAgIG52W2luZGV4T2ZOZWcxXSA9IG5cbiAgICAgICAgICBuZWlnaGJvci5ib3VuZGFyeSA9IGZhbHNlXG4gICAgICAgICAgaW50ZXJpb3IucHVzaChuZWlnaGJvcilcbiAgICAgICAgICB0b3Zpc2l0LnB1c2gobmVpZ2hib3IpXG4gICAgICAgICAgbmVpZ2hib3IubGFzdFZpc2l0ZWQgPSBuXG4gICAgICAgICAgY29udGludWVcbiAgICAgICAgfSBlbHNlIHtcbiAgICAgICAgICBuZWlnaGJvci5sYXN0VmlzaXRlZCA9IC1uXG4gICAgICAgIH1cbiAgICAgIH1cblxuICAgICAgdmFyIG5hID0gbmVpZ2hib3IuYWRqYWNlbnRcblxuICAgICAgLy9PdGhlcndpc2UsIHJlcGxhY2UgbmVpZ2hib3Igd2l0aCBuZXcgZmFjZVxuICAgICAgdmFyIHZ2ZXJ0cyA9IGNlbGxWZXJ0cy5zbGljZSgpXG4gICAgICB2YXIgdmFkaiA9IGNlbGxBZGouc2xpY2UoKVxuICAgICAgdmFyIG5jZWxsID0gbmV3IFNpbXBsZXgodnZlcnRzLCB2YWRqLCB0cnVlKVxuICAgICAgc2ltcGxpY2VzLnB1c2gobmNlbGwpXG5cbiAgICAgIC8vQ29ubmVjdCB0byBuZWlnaGJvclxuICAgICAgdmFyIG9wcG9zaXRlID0gbmEuaW5kZXhPZihjZWxsKVxuICAgICAgaWYob3Bwb3NpdGUgPCAwKSB7XG4gICAgICAgIGNvbnRpbnVlXG4gICAgICB9XG4gICAgICBuYVtvcHBvc2l0ZV0gPSBuY2VsbFxuICAgICAgdmFkaltpbmRleE9mTl0gPSBuZWlnaGJvclxuXG4gICAgICAvL0Nvbm5lY3QgdG8gY2VsbFxuICAgICAgdnZlcnRzW2ldID0gLTFcbiAgICAgIHZhZGpbaV0gPSBjZWxsXG4gICAgICBjZWxsQWRqW2ldID0gbmNlbGxcblxuICAgICAgLy9GbGlwIGZhY2V0XG4gICAgICBuY2VsbC5mbGlwKClcblxuICAgICAgLy9BZGQgdG8gZ2x1ZSBsaXN0XG4gICAgICBmb3IodmFyIGo9MDsgajw9ZDsgKytqKSB7XG4gICAgICAgIHZhciB1dSA9IHZ2ZXJ0c1tqXVxuICAgICAgICBpZih1dSA8IDAgfHwgdXUgPT09IG4pIHtcbiAgICAgICAgICBjb250aW51ZVxuICAgICAgICB9XG4gICAgICAgIHZhciBuZmFjZSA9IG5ldyBBcnJheShkLTEpXG4gICAgICAgIHZhciBucHRyID0gMFxuICAgICAgICBmb3IodmFyIGs9MDsgazw9ZDsgKytrKSB7XG4gICAgICAgICAgdmFyIHZ2ID0gdnZlcnRzW2tdXG4gICAgICAgICAgaWYodnYgPCAwIHx8IGsgPT09IGopIHtcbiAgICAgICAgICAgIGNvbnRpbnVlXG4gICAgICAgICAgfVxuICAgICAgICAgIG5mYWNlW25wdHIrK10gPSB2dlxuICAgICAgICB9XG4gICAgICAgIGdsdWVGYWNldHMucHVzaChuZXcgR2x1ZUZhY2V0KG5mYWNlLCBuY2VsbCwgaikpXG4gICAgICB9XG4gICAgfVxuICB9XG5cbiAgLy9HbHVlIGJvdW5kYXJ5IGZhY2V0cyB0b2dldGhlclxuICBnbHVlRmFjZXRzLnNvcnQoY29tcGFyZUdsdWUpXG5cbiAgZm9yKHZhciBpPTA7IGkrMTxnbHVlRmFjZXRzLmxlbmd0aDsgaSs9Mikge1xuICAgIHZhciBhID0gZ2x1ZUZhY2V0c1tpXVxuICAgIHZhciBiID0gZ2x1ZUZhY2V0c1tpKzFdXG4gICAgdmFyIGFpID0gYS5pbmRleFxuICAgIHZhciBiaSA9IGIuaW5kZXhcbiAgICBpZihhaSA8IDAgfHwgYmkgPCAwKSB7XG4gICAgICBjb250aW51ZVxuICAgIH1cbiAgICBhLmNlbGwuYWRqYWNlbnRbYS5pbmRleF0gPSBiLmNlbGxcbiAgICBiLmNlbGwuYWRqYWNlbnRbYi5pbmRleF0gPSBhLmNlbGxcbiAgfVxufVxuXG5wcm90by5pbnNlcnQgPSBmdW5jdGlvbihwb2ludCwgcmFuZG9tKSB7XG4gIC8vQWRkIHBvaW50XG4gIHZhciB2ZXJ0cyA9IHRoaXMudmVydGljZXNcbiAgdmVydHMucHVzaChwb2ludClcblxuICB2YXIgY2VsbCA9IHRoaXMud2Fsayhwb2ludCwgcmFuZG9tKVxuICBpZighY2VsbCkge1xuICAgIHJldHVyblxuICB9XG5cbiAgLy9BbGlhcyBsb2NhbCBwcm9wZXJ0aWVzXG4gIHZhciBkID0gdGhpcy5kaW1lbnNpb25cbiAgdmFyIHR1cGxlID0gdGhpcy50dXBsZVxuXG4gIC8vRGVnZW5lcmF0ZSBjYXNlOiBJZiBwb2ludCBpcyBjb3BsYW5hciB0byBjZWxsLCB0aGVuIHdhbGsgdW50aWwgd2UgZmluZCBhIG5vbi1kZWdlbmVyYXRlIGJvdW5kYXJ5XG4gIGZvcih2YXIgaT0wOyBpPD1kOyArK2kpIHtcbiAgICB2YXIgdnYgPSBjZWxsLnZlcnRpY2VzW2ldXG4gICAgaWYodnYgPCAwKSB7XG4gICAgICB0dXBsZVtpXSA9IHBvaW50XG4gICAgfSBlbHNlIHtcbiAgICAgIHR1cGxlW2ldID0gdmVydHNbdnZdXG4gICAgfVxuICB9XG4gIHZhciBvID0gdGhpcy5vcmllbnQodHVwbGUpXG4gIGlmKG8gPCAwKSB7XG4gICAgcmV0dXJuXG4gIH0gZWxzZSBpZihvID09PSAwKSB7XG4gICAgY2VsbCA9IHRoaXMuaGFuZGxlQm91bmRhcnlEZWdlbmVyYWN5KGNlbGwsIHBvaW50KVxuICAgIGlmKCFjZWxsKSB7XG4gICAgICByZXR1cm5cbiAgICB9XG4gIH1cblxuICAvL0FkZCBwZWFrc1xuICB0aGlzLmFkZFBlYWtzKHBvaW50LCBjZWxsKVxufVxuXG4vL0V4dHJhY3QgYWxsIGJvdW5kYXJ5IGNlbGxzXG5wcm90by5ib3VuZGFyeSA9IGZ1bmN0aW9uKCkge1xuICB2YXIgZCA9IHRoaXMuZGltZW5zaW9uXG4gIHZhciBib3VuZGFyeSA9IFtdXG4gIHZhciBjZWxscyA9IHRoaXMuc2ltcGxpY2VzXG4gIHZhciBuYyA9IGNlbGxzLmxlbmd0aFxuICBmb3IodmFyIGk9MDsgaTxuYzsgKytpKSB7XG4gICAgdmFyIGMgPSBjZWxsc1tpXVxuICAgIGlmKGMuYm91bmRhcnkpIHtcbiAgICAgIHZhciBiY2VsbCA9IG5ldyBBcnJheShkKVxuICAgICAgdmFyIGN2ID0gYy52ZXJ0aWNlc1xuICAgICAgdmFyIHB0ciA9IDBcbiAgICAgIHZhciBwYXJpdHkgPSAwXG4gICAgICBmb3IodmFyIGo9MDsgajw9ZDsgKytqKSB7XG4gICAgICAgIGlmKGN2W2pdID49IDApIHtcbiAgICAgICAgICBiY2VsbFtwdHIrK10gPSBjdltqXVxuICAgICAgICB9IGVsc2Uge1xuICAgICAgICAgIHBhcml0eSA9IGomMVxuICAgICAgICB9XG4gICAgICB9XG4gICAgICBpZihwYXJpdHkgPT09IChkJjEpKSB7XG4gICAgICAgIHZhciB0ID0gYmNlbGxbMF1cbiAgICAgICAgYmNlbGxbMF0gPSBiY2VsbFsxXVxuICAgICAgICBiY2VsbFsxXSA9IHRcbiAgICAgIH1cbiAgICAgIGJvdW5kYXJ5LnB1c2goYmNlbGwpXG4gICAgfVxuICB9XG4gIHJldHVybiBib3VuZGFyeVxufVxuXG5mdW5jdGlvbiBpbmNyZW1lbnRhbENvbnZleEh1bGwocG9pbnRzLCByYW5kb21TZWFyY2gpIHtcbiAgdmFyIG4gPSBwb2ludHMubGVuZ3RoXG4gIGlmKG4gPT09IDApIHtcbiAgICB0aHJvdyBuZXcgRXJyb3IoXCJNdXN0IGhhdmUgYXQgbGVhc3QgZCsxIHBvaW50c1wiKVxuICB9XG4gIHZhciBkID0gcG9pbnRzWzBdLmxlbmd0aFxuICBpZihuIDw9IGQpIHtcbiAgICB0aHJvdyBuZXcgRXJyb3IoXCJNdXN0IGlucHV0IGF0IGxlYXN0IGQrMSBwb2ludHNcIilcbiAgfVxuXG4gIC8vRklYTUU6IFRoaXMgY291bGQgYmUgZGVnZW5lcmF0ZSwgYnV0IG5lZWQgdG8gc2VsZWN0IGQrMSBub24tY29wbGFuYXIgcG9pbnRzIHRvIGJvb3RzdHJhcCBwcm9jZXNzXG4gIHZhciBpbml0aWFsU2ltcGxleCA9IHBvaW50cy5zbGljZSgwLCBkKzEpXG5cbiAgLy9NYWtlIHN1cmUgaW5pdGlhbCBzaW1wbGV4IGlzIHBvc2l0aXZlbHkgb3JpZW50ZWRcbiAgdmFyIG8gPSBvcmllbnQuYXBwbHkodm9pZCAwLCBpbml0aWFsU2ltcGxleClcbiAgaWYobyA9PT0gMCkge1xuICAgIHRocm93IG5ldyBFcnJvcihcIklucHV0IG5vdCBpbiBnZW5lcmFsIHBvc2l0aW9uXCIpXG4gIH1cbiAgdmFyIGluaXRpYWxDb29yZHMgPSBuZXcgQXJyYXkoZCsxKVxuICBmb3IodmFyIGk9MDsgaTw9ZDsgKytpKSB7XG4gICAgaW5pdGlhbENvb3Jkc1tpXSA9IGlcbiAgfVxuICBpZihvIDwgMCkge1xuICAgIGluaXRpYWxDb29yZHNbMF0gPSAxXG4gICAgaW5pdGlhbENvb3Jkc1sxXSA9IDBcbiAgfVxuXG4gIC8vQ3JlYXRlIGluaXRpYWwgdG9wb2xvZ2ljYWwgaW5kZXgsIGdsdWUgcG9pbnRlcnMgdG9nZXRoZXIgKGtpbmQgb2YgbWVzc3kpXG4gIHZhciBpbml0aWFsQ2VsbCA9IG5ldyBTaW1wbGV4KGluaXRpYWxDb29yZHMsIG5ldyBBcnJheShkKzEpLCBmYWxzZSlcbiAgdmFyIGJvdW5kYXJ5ID0gaW5pdGlhbENlbGwuYWRqYWNlbnRcbiAgdmFyIGxpc3QgPSBuZXcgQXJyYXkoZCsyKVxuICBmb3IodmFyIGk9MDsgaTw9ZDsgKytpKSB7XG4gICAgdmFyIHZlcnRzID0gaW5pdGlhbENvb3Jkcy5zbGljZSgpXG4gICAgZm9yKHZhciBqPTA7IGo8PWQ7ICsraikge1xuICAgICAgaWYoaiA9PT0gaSkge1xuICAgICAgICB2ZXJ0c1tqXSA9IC0xXG4gICAgICB9XG4gICAgfVxuICAgIHZhciB0ID0gdmVydHNbMF1cbiAgICB2ZXJ0c1swXSA9IHZlcnRzWzFdXG4gICAgdmVydHNbMV0gPSB0XG4gICAgdmFyIGNlbGwgPSBuZXcgU2ltcGxleCh2ZXJ0cywgbmV3IEFycmF5KGQrMSksIHRydWUpXG4gICAgYm91bmRhcnlbaV0gPSBjZWxsXG4gICAgbGlzdFtpXSA9IGNlbGxcbiAgfVxuICBsaXN0W2QrMV0gPSBpbml0aWFsQ2VsbFxuICBmb3IodmFyIGk9MDsgaTw9ZDsgKytpKSB7XG4gICAgdmFyIHZlcnRzID0gYm91bmRhcnlbaV0udmVydGljZXNcbiAgICB2YXIgYWRqID0gYm91bmRhcnlbaV0uYWRqYWNlbnRcbiAgICBmb3IodmFyIGo9MDsgajw9ZDsgKytqKSB7XG4gICAgICB2YXIgdiA9IHZlcnRzW2pdXG4gICAgICBpZih2IDwgMCkge1xuICAgICAgICBhZGpbal0gPSBpbml0aWFsQ2VsbFxuICAgICAgICBjb250aW51ZVxuICAgICAgfVxuICAgICAgZm9yKHZhciBrPTA7IGs8PWQ7ICsraykge1xuICAgICAgICBpZihib3VuZGFyeVtrXS52ZXJ0aWNlcy5pbmRleE9mKHYpIDwgMCkge1xuICAgICAgICAgIGFkaltqXSA9IGJvdW5kYXJ5W2tdXG4gICAgICAgIH1cbiAgICAgIH1cbiAgICB9XG4gIH1cblxuICAvL0luaXRpYWxpemUgdHJpYW5nbGVzXG4gIHZhciB0cmlhbmdsZXMgPSBuZXcgVHJpYW5ndWxhdGlvbihkLCBpbml0aWFsU2ltcGxleCwgbGlzdClcblxuICAvL0luc2VydCByZW1haW5pbmcgcG9pbnRzXG4gIHZhciB1c2VSYW5kb20gPSAhIXJhbmRvbVNlYXJjaFxuICBmb3IodmFyIGk9ZCsxOyBpPG47ICsraSkge1xuICAgIHRyaWFuZ2xlcy5pbnNlcnQocG9pbnRzW2ldLCB1c2VSYW5kb20pXG4gIH1cbiAgXG4gIC8vRXh0cmFjdCBib3VuZGFyeSBjZWxsc1xuICByZXR1cm4gdHJpYW5nbGVzLmJvdW5kYXJ5KClcbn0iLCJcInVzZSBzdHJpY3RcIlxuXG5tb2R1bGUuZXhwb3J0cyA9IGZhc3RUd29TdW1cblxuZnVuY3Rpb24gZmFzdFR3b1N1bShhLCBiLCByZXN1bHQpIHtcblx0dmFyIHggPSBhICsgYlxuXHR2YXIgYnYgPSB4IC0gYVxuXHR2YXIgYXYgPSB4IC0gYnZcblx0dmFyIGJyID0gYiAtIGJ2XG5cdHZhciBhciA9IGEgLSBhdlxuXHRpZihyZXN1bHQpIHtcblx0XHRyZXN1bHRbMF0gPSBhciArIGJyXG5cdFx0cmVzdWx0WzFdID0geFxuXHRcdHJldHVybiByZXN1bHRcblx0fVxuXHRyZXR1cm4gW2FyK2JyLCB4XVxufSIsIlwidXNlIHN0cmljdFwiXG5cbnZhciB0d29Qcm9kdWN0ID0gcmVxdWlyZShcInR3by1wcm9kdWN0XCIpXG52YXIgdHdvU3VtID0gcmVxdWlyZShcInR3by1zdW1cIilcblxubW9kdWxlLmV4cG9ydHMgPSBzY2FsZUxpbmVhckV4cGFuc2lvblxuXG5mdW5jdGlvbiBzY2FsZUxpbmVhckV4cGFuc2lvbihlLCBzY2FsZSkge1xuICB2YXIgbiA9IGUubGVuZ3RoXG4gIGlmKG4gPT09IDEpIHtcbiAgICB2YXIgdHMgPSB0d29Qcm9kdWN0KGVbMF0sIHNjYWxlKVxuICAgIGlmKHRzWzBdKSB7XG4gICAgICByZXR1cm4gdHNcbiAgICB9XG4gICAgcmV0dXJuIFsgdHNbMV0gXVxuICB9XG4gIHZhciBnID0gbmV3IEFycmF5KDIgKiBuKVxuICB2YXIgcSA9IFswLjEsIDAuMV1cbiAgdmFyIHQgPSBbMC4xLCAwLjFdXG4gIHZhciBjb3VudCA9IDBcbiAgdHdvUHJvZHVjdChlWzBdLCBzY2FsZSwgcSlcbiAgaWYocVswXSkge1xuICAgIGdbY291bnQrK10gPSBxWzBdXG4gIH1cbiAgZm9yKHZhciBpPTE7IGk8bjsgKytpKSB7XG4gICAgdHdvUHJvZHVjdChlW2ldLCBzY2FsZSwgdClcbiAgICB2YXIgcHEgPSBxWzFdXG4gICAgdHdvU3VtKHBxLCB0WzBdLCBxKVxuICAgIGlmKHFbMF0pIHtcbiAgICAgIGdbY291bnQrK10gPSBxWzBdXG4gICAgfVxuICAgIHZhciBhID0gdFsxXVxuICAgIHZhciBiID0gcVsxXVxuICAgIHZhciB4ID0gYSArIGJcbiAgICB2YXIgYnYgPSB4IC0gYVxuICAgIHZhciB5ID0gYiAtIGJ2XG4gICAgcVsxXSA9IHhcbiAgICBpZih5KSB7XG4gICAgICBnW2NvdW50KytdID0geVxuICAgIH1cbiAgfVxuICBpZihxWzFdKSB7XG4gICAgZ1tjb3VudCsrXSA9IHFbMV1cbiAgfVxuICBpZihjb3VudCA9PT0gMCkge1xuICAgIGdbY291bnQrK10gPSAwLjBcbiAgfVxuICBnLmxlbmd0aCA9IGNvdW50XG4gIHJldHVybiBnXG59IiwiXCJ1c2Ugc3RyaWN0XCJcblxubW9kdWxlLmV4cG9ydHMgPSByb2J1c3RTdWJ0cmFjdFxuXG4vL0Vhc3kgY2FzZTogQWRkIHR3byBzY2FsYXJzXG5mdW5jdGlvbiBzY2FsYXJTY2FsYXIoYSwgYikge1xuICB2YXIgeCA9IGEgKyBiXG4gIHZhciBidiA9IHggLSBhXG4gIHZhciBhdiA9IHggLSBidlxuICB2YXIgYnIgPSBiIC0gYnZcbiAgdmFyIGFyID0gYSAtIGF2XG4gIHZhciB5ID0gYXIgKyBiclxuICBpZih5KSB7XG4gICAgcmV0dXJuIFt5LCB4XVxuICB9XG4gIHJldHVybiBbeF1cbn1cblxuZnVuY3Rpb24gcm9idXN0U3VidHJhY3QoZSwgZikge1xuICB2YXIgbmUgPSBlLmxlbmd0aHwwXG4gIHZhciBuZiA9IGYubGVuZ3RofDBcbiAgaWYobmUgPT09IDEgJiYgbmYgPT09IDEpIHtcbiAgICByZXR1cm4gc2NhbGFyU2NhbGFyKGVbMF0sIC1mWzBdKVxuICB9XG4gIHZhciBuID0gbmUgKyBuZlxuICB2YXIgZyA9IG5ldyBBcnJheShuKVxuICB2YXIgY291bnQgPSAwXG4gIHZhciBlcHRyID0gMFxuICB2YXIgZnB0ciA9IDBcbiAgdmFyIGFicyA9IE1hdGguYWJzXG4gIHZhciBlaSA9IGVbZXB0cl1cbiAgdmFyIGVhID0gYWJzKGVpKVxuICB2YXIgZmkgPSAtZltmcHRyXVxuICB2YXIgZmEgPSBhYnMoZmkpXG4gIHZhciBhLCBiXG4gIGlmKGVhIDwgZmEpIHtcbiAgICBiID0gZWlcbiAgICBlcHRyICs9IDFcbiAgICBpZihlcHRyIDwgbmUpIHtcbiAgICAgIGVpID0gZVtlcHRyXVxuICAgICAgZWEgPSBhYnMoZWkpXG4gICAgfVxuICB9IGVsc2Uge1xuICAgIGIgPSBmaVxuICAgIGZwdHIgKz0gMVxuICAgIGlmKGZwdHIgPCBuZikge1xuICAgICAgZmkgPSAtZltmcHRyXVxuICAgICAgZmEgPSBhYnMoZmkpXG4gICAgfVxuICB9XG4gIGlmKChlcHRyIDwgbmUgJiYgZWEgPCBmYSkgfHwgKGZwdHIgPj0gbmYpKSB7XG4gICAgYSA9IGVpXG4gICAgZXB0ciArPSAxXG4gICAgaWYoZXB0ciA8IG5lKSB7XG4gICAgICBlaSA9IGVbZXB0cl1cbiAgICAgIGVhID0gYWJzKGVpKVxuICAgIH1cbiAgfSBlbHNlIHtcbiAgICBhID0gZmlcbiAgICBmcHRyICs9IDFcbiAgICBpZihmcHRyIDwgbmYpIHtcbiAgICAgIGZpID0gLWZbZnB0cl1cbiAgICAgIGZhID0gYWJzKGZpKVxuICAgIH1cbiAgfVxuICB2YXIgeCA9IGEgKyBiXG4gIHZhciBidiA9IHggLSBhXG4gIHZhciB5ID0gYiAtIGJ2XG4gIHZhciBxMCA9IHlcbiAgdmFyIHExID0geFxuICB2YXIgX3gsIF9idiwgX2F2LCBfYnIsIF9hclxuICB3aGlsZShlcHRyIDwgbmUgJiYgZnB0ciA8IG5mKSB7XG4gICAgaWYoZWEgPCBmYSkge1xuICAgICAgYSA9IGVpXG4gICAgICBlcHRyICs9IDFcbiAgICAgIGlmKGVwdHIgPCBuZSkge1xuICAgICAgICBlaSA9IGVbZXB0cl1cbiAgICAgICAgZWEgPSBhYnMoZWkpXG4gICAgICB9XG4gICAgfSBlbHNlIHtcbiAgICAgIGEgPSBmaVxuICAgICAgZnB0ciArPSAxXG4gICAgICBpZihmcHRyIDwgbmYpIHtcbiAgICAgICAgZmkgPSAtZltmcHRyXVxuICAgICAgICBmYSA9IGFicyhmaSlcbiAgICAgIH1cbiAgICB9XG4gICAgYiA9IHEwXG4gICAgeCA9IGEgKyBiXG4gICAgYnYgPSB4IC0gYVxuICAgIHkgPSBiIC0gYnZcbiAgICBpZih5KSB7XG4gICAgICBnW2NvdW50KytdID0geVxuICAgIH1cbiAgICBfeCA9IHExICsgeFxuICAgIF9idiA9IF94IC0gcTFcbiAgICBfYXYgPSBfeCAtIF9idlxuICAgIF9iciA9IHggLSBfYnZcbiAgICBfYXIgPSBxMSAtIF9hdlxuICAgIHEwID0gX2FyICsgX2JyXG4gICAgcTEgPSBfeFxuICB9XG4gIHdoaWxlKGVwdHIgPCBuZSkge1xuICAgIGEgPSBlaVxuICAgIGIgPSBxMFxuICAgIHggPSBhICsgYlxuICAgIGJ2ID0geCAtIGFcbiAgICB5ID0gYiAtIGJ2XG4gICAgaWYoeSkge1xuICAgICAgZ1tjb3VudCsrXSA9IHlcbiAgICB9XG4gICAgX3ggPSBxMSArIHhcbiAgICBfYnYgPSBfeCAtIHExXG4gICAgX2F2ID0gX3ggLSBfYnZcbiAgICBfYnIgPSB4IC0gX2J2XG4gICAgX2FyID0gcTEgLSBfYXZcbiAgICBxMCA9IF9hciArIF9iclxuICAgIHExID0gX3hcbiAgICBlcHRyICs9IDFcbiAgICBpZihlcHRyIDwgbmUpIHtcbiAgICAgIGVpID0gZVtlcHRyXVxuICAgIH1cbiAgfVxuICB3aGlsZShmcHRyIDwgbmYpIHtcbiAgICBhID0gZmlcbiAgICBiID0gcTBcbiAgICB4ID0gYSArIGJcbiAgICBidiA9IHggLSBhXG4gICAgeSA9IGIgLSBidlxuICAgIGlmKHkpIHtcbiAgICAgIGdbY291bnQrK10gPSB5XG4gICAgfSBcbiAgICBfeCA9IHExICsgeFxuICAgIF9idiA9IF94IC0gcTFcbiAgICBfYXYgPSBfeCAtIF9idlxuICAgIF9iciA9IHggLSBfYnZcbiAgICBfYXIgPSBxMSAtIF9hdlxuICAgIHEwID0gX2FyICsgX2JyXG4gICAgcTEgPSBfeFxuICAgIGZwdHIgKz0gMVxuICAgIGlmKGZwdHIgPCBuZikge1xuICAgICAgZmkgPSAtZltmcHRyXVxuICAgIH1cbiAgfVxuICBpZihxMCkge1xuICAgIGdbY291bnQrK10gPSBxMFxuICB9XG4gIGlmKHExKSB7XG4gICAgZ1tjb3VudCsrXSA9IHExXG4gIH1cbiAgaWYoIWNvdW50KSB7XG4gICAgZ1tjb3VudCsrXSA9IDAuMCAgXG4gIH1cbiAgZy5sZW5ndGggPSBjb3VudFxuICByZXR1cm4gZ1xufSIsIlwidXNlIHN0cmljdFwiXG5cbm1vZHVsZS5leHBvcnRzID0gbGluZWFyRXhwYW5zaW9uU3VtXG5cbi8vRWFzeSBjYXNlOiBBZGQgdHdvIHNjYWxhcnNcbmZ1bmN0aW9uIHNjYWxhclNjYWxhcihhLCBiKSB7XG4gIHZhciB4ID0gYSArIGJcbiAgdmFyIGJ2ID0geCAtIGFcbiAgdmFyIGF2ID0geCAtIGJ2XG4gIHZhciBiciA9IGIgLSBidlxuICB2YXIgYXIgPSBhIC0gYXZcbiAgdmFyIHkgPSBhciArIGJyXG4gIGlmKHkpIHtcbiAgICByZXR1cm4gW3ksIHhdXG4gIH1cbiAgcmV0dXJuIFt4XVxufVxuXG5mdW5jdGlvbiBsaW5lYXJFeHBhbnNpb25TdW0oZSwgZikge1xuICB2YXIgbmUgPSBlLmxlbmd0aHwwXG4gIHZhciBuZiA9IGYubGVuZ3RofDBcbiAgaWYobmUgPT09IDEgJiYgbmYgPT09IDEpIHtcbiAgICByZXR1cm4gc2NhbGFyU2NhbGFyKGVbMF0sIGZbMF0pXG4gIH1cbiAgdmFyIG4gPSBuZSArIG5mXG4gIHZhciBnID0gbmV3IEFycmF5KG4pXG4gIHZhciBjb3VudCA9IDBcbiAgdmFyIGVwdHIgPSAwXG4gIHZhciBmcHRyID0gMFxuICB2YXIgYWJzID0gTWF0aC5hYnNcbiAgdmFyIGVpID0gZVtlcHRyXVxuICB2YXIgZWEgPSBhYnMoZWkpXG4gIHZhciBmaSA9IGZbZnB0cl1cbiAgdmFyIGZhID0gYWJzKGZpKVxuICB2YXIgYSwgYlxuICBpZihlYSA8IGZhKSB7XG4gICAgYiA9IGVpXG4gICAgZXB0ciArPSAxXG4gICAgaWYoZXB0ciA8IG5lKSB7XG4gICAgICBlaSA9IGVbZXB0cl1cbiAgICAgIGVhID0gYWJzKGVpKVxuICAgIH1cbiAgfSBlbHNlIHtcbiAgICBiID0gZmlcbiAgICBmcHRyICs9IDFcbiAgICBpZihmcHRyIDwgbmYpIHtcbiAgICAgIGZpID0gZltmcHRyXVxuICAgICAgZmEgPSBhYnMoZmkpXG4gICAgfVxuICB9XG4gIGlmKChlcHRyIDwgbmUgJiYgZWEgPCBmYSkgfHwgKGZwdHIgPj0gbmYpKSB7XG4gICAgYSA9IGVpXG4gICAgZXB0ciArPSAxXG4gICAgaWYoZXB0ciA8IG5lKSB7XG4gICAgICBlaSA9IGVbZXB0cl1cbiAgICAgIGVhID0gYWJzKGVpKVxuICAgIH1cbiAgfSBlbHNlIHtcbiAgICBhID0gZmlcbiAgICBmcHRyICs9IDFcbiAgICBpZihmcHRyIDwgbmYpIHtcbiAgICAgIGZpID0gZltmcHRyXVxuICAgICAgZmEgPSBhYnMoZmkpXG4gICAgfVxuICB9XG4gIHZhciB4ID0gYSArIGJcbiAgdmFyIGJ2ID0geCAtIGFcbiAgdmFyIHkgPSBiIC0gYnZcbiAgdmFyIHEwID0geVxuICB2YXIgcTEgPSB4XG4gIHZhciBfeCwgX2J2LCBfYXYsIF9iciwgX2FyXG4gIHdoaWxlKGVwdHIgPCBuZSAmJiBmcHRyIDwgbmYpIHtcbiAgICBpZihlYSA8IGZhKSB7XG4gICAgICBhID0gZWlcbiAgICAgIGVwdHIgKz0gMVxuICAgICAgaWYoZXB0ciA8IG5lKSB7XG4gICAgICAgIGVpID0gZVtlcHRyXVxuICAgICAgICBlYSA9IGFicyhlaSlcbiAgICAgIH1cbiAgICB9IGVsc2Uge1xuICAgICAgYSA9IGZpXG4gICAgICBmcHRyICs9IDFcbiAgICAgIGlmKGZwdHIgPCBuZikge1xuICAgICAgICBmaSA9IGZbZnB0cl1cbiAgICAgICAgZmEgPSBhYnMoZmkpXG4gICAgICB9XG4gICAgfVxuICAgIGIgPSBxMFxuICAgIHggPSBhICsgYlxuICAgIGJ2ID0geCAtIGFcbiAgICB5ID0gYiAtIGJ2XG4gICAgaWYoeSkge1xuICAgICAgZ1tjb3VudCsrXSA9IHlcbiAgICB9XG4gICAgX3ggPSBxMSArIHhcbiAgICBfYnYgPSBfeCAtIHExXG4gICAgX2F2ID0gX3ggLSBfYnZcbiAgICBfYnIgPSB4IC0gX2J2XG4gICAgX2FyID0gcTEgLSBfYXZcbiAgICBxMCA9IF9hciArIF9iclxuICAgIHExID0gX3hcbiAgfVxuICB3aGlsZShlcHRyIDwgbmUpIHtcbiAgICBhID0gZWlcbiAgICBiID0gcTBcbiAgICB4ID0gYSArIGJcbiAgICBidiA9IHggLSBhXG4gICAgeSA9IGIgLSBidlxuICAgIGlmKHkpIHtcbiAgICAgIGdbY291bnQrK10gPSB5XG4gICAgfVxuICAgIF94ID0gcTEgKyB4XG4gICAgX2J2ID0gX3ggLSBxMVxuICAgIF9hdiA9IF94IC0gX2J2XG4gICAgX2JyID0geCAtIF9idlxuICAgIF9hciA9IHExIC0gX2F2XG4gICAgcTAgPSBfYXIgKyBfYnJcbiAgICBxMSA9IF94XG4gICAgZXB0ciArPSAxXG4gICAgaWYoZXB0ciA8IG5lKSB7XG4gICAgICBlaSA9IGVbZXB0cl1cbiAgICB9XG4gIH1cbiAgd2hpbGUoZnB0ciA8IG5mKSB7XG4gICAgYSA9IGZpXG4gICAgYiA9IHEwXG4gICAgeCA9IGEgKyBiXG4gICAgYnYgPSB4IC0gYVxuICAgIHkgPSBiIC0gYnZcbiAgICBpZih5KSB7XG4gICAgICBnW2NvdW50KytdID0geVxuICAgIH0gXG4gICAgX3ggPSBxMSArIHhcbiAgICBfYnYgPSBfeCAtIHExXG4gICAgX2F2ID0gX3ggLSBfYnZcbiAgICBfYnIgPSB4IC0gX2J2XG4gICAgX2FyID0gcTEgLSBfYXZcbiAgICBxMCA9IF9hciArIF9iclxuICAgIHExID0gX3hcbiAgICBmcHRyICs9IDFcbiAgICBpZihmcHRyIDwgbmYpIHtcbiAgICAgIGZpID0gZltmcHRyXVxuICAgIH1cbiAgfVxuICBpZihxMCkge1xuICAgIGdbY291bnQrK10gPSBxMFxuICB9XG4gIGlmKHExKSB7XG4gICAgZ1tjb3VudCsrXSA9IHExXG4gIH1cbiAgaWYoIWNvdW50KSB7XG4gICAgZ1tjb3VudCsrXSA9IDAuMCAgXG4gIH1cbiAgZy5sZW5ndGggPSBjb3VudFxuICByZXR1cm4gZ1xufSIsIlwidXNlIHN0cmljdFwiXG5cbm1vZHVsZS5leHBvcnRzID0gdHdvUHJvZHVjdFxuXG52YXIgU1BMSVRURVIgPSArKE1hdGgucG93KDIsIDI3KSArIDEuMClcblxuZnVuY3Rpb24gdHdvUHJvZHVjdChhLCBiLCByZXN1bHQpIHtcbiAgdmFyIHggPSBhICogYlxuXG4gIHZhciBjID0gU1BMSVRURVIgKiBhXG4gIHZhciBhYmlnID0gYyAtIGFcbiAgdmFyIGFoaSA9IGMgLSBhYmlnXG4gIHZhciBhbG8gPSBhIC0gYWhpXG5cbiAgdmFyIGQgPSBTUExJVFRFUiAqIGJcbiAgdmFyIGJiaWcgPSBkIC0gYlxuICB2YXIgYmhpID0gZCAtIGJiaWdcbiAgdmFyIGJsbyA9IGIgLSBiaGlcblxuICB2YXIgZXJyMSA9IHggLSAoYWhpICogYmhpKVxuICB2YXIgZXJyMiA9IGVycjEgLSAoYWxvICogYmhpKVxuICB2YXIgZXJyMyA9IGVycjIgLSAoYWhpICogYmxvKVxuXG4gIHZhciB5ID0gYWxvICogYmxvIC0gZXJyM1xuXG4gIGlmKHJlc3VsdCkge1xuICAgIHJlc3VsdFswXSA9IHlcbiAgICByZXN1bHRbMV0gPSB4XG4gICAgcmV0dXJuIHJlc3VsdFxuICB9XG5cbiAgcmV0dXJuIFsgeSwgeCBdXG59IiwiXCJ1c2Ugc3RyaWN0XCJcblxudmFyIHR3b1Byb2R1Y3QgPSByZXF1aXJlKFwidHdvLXByb2R1Y3RcIilcbnZhciByb2J1c3RTdW0gPSByZXF1aXJlKFwicm9idXN0LXN1bVwiKVxudmFyIHJvYnVzdFNjYWxlID0gcmVxdWlyZShcInJvYnVzdC1zY2FsZVwiKVxudmFyIHJvYnVzdFN1YnRyYWN0ID0gcmVxdWlyZShcInJvYnVzdC1zdWJ0cmFjdFwiKVxuXG52YXIgTlVNX0VYUEFORCA9IDVcblxudmFyIEVQU0lMT04gICAgID0gMS4xMTAyMjMwMjQ2MjUxNTY1ZS0xNlxudmFyIEVSUkJPVU5EMyAgID0gKDMuMCArIDE2LjAgKiBFUFNJTE9OKSAqIEVQU0lMT05cbnZhciBFUlJCT1VORDQgICA9ICg3LjAgKyA1Ni4wICogRVBTSUxPTikgKiBFUFNJTE9OXG5cbmZ1bmN0aW9uIGNvZmFjdG9yKG0sIGMpIHtcbiAgdmFyIHJlc3VsdCA9IG5ldyBBcnJheShtLmxlbmd0aC0xKVxuICBmb3IodmFyIGk9MTsgaTxtLmxlbmd0aDsgKytpKSB7XG4gICAgdmFyIHIgPSByZXN1bHRbaS0xXSA9IG5ldyBBcnJheShtLmxlbmd0aC0xKVxuICAgIGZvcih2YXIgaj0wLGs9MDsgajxtLmxlbmd0aDsgKytqKSB7XG4gICAgICBpZihqID09PSBjKSB7XG4gICAgICAgIGNvbnRpbnVlXG4gICAgICB9XG4gICAgICByW2srK10gPSBtW2ldW2pdXG4gICAgfVxuICB9XG4gIHJldHVybiByZXN1bHRcbn1cblxuZnVuY3Rpb24gbWF0cml4KG4pIHtcbiAgdmFyIHJlc3VsdCA9IG5ldyBBcnJheShuKVxuICBmb3IodmFyIGk9MDsgaTxuOyArK2kpIHtcbiAgICByZXN1bHRbaV0gPSBuZXcgQXJyYXkobilcbiAgICBmb3IodmFyIGo9MDsgajxuOyArK2opIHtcbiAgICAgIHJlc3VsdFtpXVtqXSA9IFtcIm1cIiwgaiwgXCJbXCIsIChuLWktMSksIFwiXVwiXS5qb2luKFwiXCIpXG4gICAgfVxuICB9XG4gIHJldHVybiByZXN1bHRcbn1cblxuZnVuY3Rpb24gc2lnbihuKSB7XG4gIGlmKG4gJiAxKSB7XG4gICAgcmV0dXJuIFwiLVwiXG4gIH1cbiAgcmV0dXJuIFwiXCJcbn1cblxuZnVuY3Rpb24gZ2VuZXJhdGVTdW0oZXhwcikge1xuICBpZihleHByLmxlbmd0aCA9PT0gMSkge1xuICAgIHJldHVybiBleHByWzBdXG4gIH0gZWxzZSBpZihleHByLmxlbmd0aCA9PT0gMikge1xuICAgIHJldHVybiBbXCJzdW0oXCIsIGV4cHJbMF0sIFwiLFwiLCBleHByWzFdLCBcIilcIl0uam9pbihcIlwiKVxuICB9IGVsc2Uge1xuICAgIHZhciBtID0gZXhwci5sZW5ndGg+PjFcbiAgICByZXR1cm4gW1wic3VtKFwiLCBnZW5lcmF0ZVN1bShleHByLnNsaWNlKDAsIG0pKSwgXCIsXCIsIGdlbmVyYXRlU3VtKGV4cHIuc2xpY2UobSkpLCBcIilcIl0uam9pbihcIlwiKVxuICB9XG59XG5cbmZ1bmN0aW9uIGRldGVybWluYW50KG0pIHtcbiAgaWYobS5sZW5ndGggPT09IDIpIHtcbiAgICByZXR1cm4gW1tcInN1bShwcm9kKFwiLCBtWzBdWzBdLCBcIixcIiwgbVsxXVsxXSwgXCIpLHByb2QoLVwiLCBtWzBdWzFdLCBcIixcIiwgbVsxXVswXSwgXCIpKVwiXS5qb2luKFwiXCIpXVxuICB9IGVsc2Uge1xuICAgIHZhciBleHByID0gW11cbiAgICBmb3IodmFyIGk9MDsgaTxtLmxlbmd0aDsgKytpKSB7XG4gICAgICBleHByLnB1c2goW1wic2NhbGUoXCIsIGdlbmVyYXRlU3VtKGRldGVybWluYW50KGNvZmFjdG9yKG0sIGkpKSksIFwiLFwiLCBzaWduKGkpLCBtWzBdW2ldLCBcIilcIl0uam9pbihcIlwiKSlcbiAgICB9XG4gICAgcmV0dXJuIGV4cHJcbiAgfVxufVxuXG5mdW5jdGlvbiBvcmllbnRhdGlvbihuKSB7XG4gIHZhciBwb3MgPSBbXVxuICB2YXIgbmVnID0gW11cbiAgdmFyIG0gPSBtYXRyaXgobilcbiAgdmFyIGFyZ3MgPSBbXVxuICBmb3IodmFyIGk9MDsgaTxuOyArK2kpIHtcbiAgICBpZigoaSYxKT09PTApIHtcbiAgICAgIHBvcy5wdXNoLmFwcGx5KHBvcywgZGV0ZXJtaW5hbnQoY29mYWN0b3IobSwgaSkpKVxuICAgIH0gZWxzZSB7XG4gICAgICBuZWcucHVzaC5hcHBseShuZWcsIGRldGVybWluYW50KGNvZmFjdG9yKG0sIGkpKSlcbiAgICB9XG4gICAgYXJncy5wdXNoKFwibVwiICsgaSlcbiAgfVxuICB2YXIgcG9zRXhwciA9IGdlbmVyYXRlU3VtKHBvcylcbiAgdmFyIG5lZ0V4cHIgPSBnZW5lcmF0ZVN1bShuZWcpXG4gIHZhciBmdW5jTmFtZSA9IFwib3JpZW50YXRpb25cIiArIG4gKyBcIkV4YWN0XCJcbiAgdmFyIGNvZGUgPSBbXCJmdW5jdGlvbiBcIiwgZnVuY05hbWUsIFwiKFwiLCBhcmdzLmpvaW4oKSwgXCIpe3ZhciBwPVwiLCBwb3NFeHByLCBcIixuPVwiLCBuZWdFeHByLCBcIixkPXN1YihwLG4pO1xcXG5yZXR1cm4gZFtkLmxlbmd0aC0xXTt9O3JldHVybiBcIiwgZnVuY05hbWVdLmpvaW4oXCJcIilcbiAgdmFyIHByb2MgPSBuZXcgRnVuY3Rpb24oXCJzdW1cIiwgXCJwcm9kXCIsIFwic2NhbGVcIiwgXCJzdWJcIiwgY29kZSlcbiAgcmV0dXJuIHByb2Mocm9idXN0U3VtLCB0d29Qcm9kdWN0LCByb2J1c3RTY2FsZSwgcm9idXN0U3VidHJhY3QpXG59XG5cbnZhciBvcmllbnRhdGlvbjNFeGFjdCA9IG9yaWVudGF0aW9uKDMpXG52YXIgb3JpZW50YXRpb240RXhhY3QgPSBvcmllbnRhdGlvbig0KVxuXG52YXIgQ0FDSEVEID0gW1xuICBmdW5jdGlvbiBvcmllbnRhdGlvbjAoKSB7IHJldHVybiAwIH0sXG4gIGZ1bmN0aW9uIG9yaWVudGF0aW9uMSgpIHsgcmV0dXJuIDAgfSxcbiAgZnVuY3Rpb24gb3JpZW50YXRpb24yKGEsIGIpIHsgXG4gICAgcmV0dXJuIGJbMF0gLSBhWzBdXG4gIH0sXG4gIGZ1bmN0aW9uIG9yaWVudGF0aW9uMyhhLCBiLCBjKSB7XG4gICAgdmFyIGwgPSAoYVsxXSAtIGNbMV0pICogKGJbMF0gLSBjWzBdKVxuICAgIHZhciByID0gKGFbMF0gLSBjWzBdKSAqIChiWzFdIC0gY1sxXSlcbiAgICB2YXIgZGV0ID0gbCAtIHJcbiAgICB2YXIgc1xuICAgIGlmKGwgPiAwKSB7XG4gICAgICBpZihyIDw9IDApIHtcbiAgICAgICAgcmV0dXJuIGRldFxuICAgICAgfSBlbHNlIHtcbiAgICAgICAgcyA9IGwgKyByXG4gICAgICB9XG4gICAgfSBlbHNlIGlmKGwgPCAwKSB7XG4gICAgICBpZihyID49IDApIHtcbiAgICAgICAgcmV0dXJuIGRldFxuICAgICAgfSBlbHNlIHtcbiAgICAgICAgcyA9IC0obCArIHIpXG4gICAgICB9XG4gICAgfSBlbHNlIHtcbiAgICAgIHJldHVybiBkZXRcbiAgICB9XG4gICAgdmFyIHRvbCA9IEVSUkJPVU5EMyAqIHNcbiAgICBpZihkZXQgPj0gdG9sIHx8IGRldCA8PSAtdG9sKSB7XG4gICAgICByZXR1cm4gZGV0XG4gICAgfVxuICAgIHJldHVybiBvcmllbnRhdGlvbjNFeGFjdChhLCBiLCBjKVxuICB9LFxuICBmdW5jdGlvbiBvcmllbnRhdGlvbjQoYSxiLGMsZCkge1xuICAgIHZhciBhZHggPSBhWzBdIC0gZFswXVxuICAgIHZhciBiZHggPSBiWzBdIC0gZFswXVxuICAgIHZhciBjZHggPSBjWzBdIC0gZFswXVxuICAgIHZhciBhZHkgPSBhWzFdIC0gZFsxXVxuICAgIHZhciBiZHkgPSBiWzFdIC0gZFsxXVxuICAgIHZhciBjZHkgPSBjWzFdIC0gZFsxXVxuICAgIHZhciBhZHogPSBhWzJdIC0gZFsyXVxuICAgIHZhciBiZHogPSBiWzJdIC0gZFsyXVxuICAgIHZhciBjZHogPSBjWzJdIC0gZFsyXVxuICAgIHZhciBiZHhjZHkgPSBiZHggKiBjZHlcbiAgICB2YXIgY2R4YmR5ID0gY2R4ICogYmR5XG4gICAgdmFyIGNkeGFkeSA9IGNkeCAqIGFkeVxuICAgIHZhciBhZHhjZHkgPSBhZHggKiBjZHlcbiAgICB2YXIgYWR4YmR5ID0gYWR4ICogYmR5XG4gICAgdmFyIGJkeGFkeSA9IGJkeCAqIGFkeVxuICAgIHZhciBkZXQgPSBhZHogKiAoYmR4Y2R5IC0gY2R4YmR5KSBcbiAgICAgICAgICAgICsgYmR6ICogKGNkeGFkeSAtIGFkeGNkeSlcbiAgICAgICAgICAgICsgY2R6ICogKGFkeGJkeSAtIGJkeGFkeSlcbiAgICB2YXIgcGVybWFuZW50ID0gKE1hdGguYWJzKGJkeGNkeSkgKyBNYXRoLmFicyhjZHhiZHkpKSAqIE1hdGguYWJzKGFkeilcbiAgICAgICAgICAgICAgICAgICsgKE1hdGguYWJzKGNkeGFkeSkgKyBNYXRoLmFicyhhZHhjZHkpKSAqIE1hdGguYWJzKGJkeilcbiAgICAgICAgICAgICAgICAgICsgKE1hdGguYWJzKGFkeGJkeSkgKyBNYXRoLmFicyhiZHhhZHkpKSAqIE1hdGguYWJzKGNkeilcbiAgICB2YXIgdG9sID0gRVJSQk9VTkQ0ICogcGVybWFuZW50XG4gICAgaWYgKChkZXQgPiB0b2wpIHx8ICgtZGV0ID4gdG9sKSkge1xuICAgICAgcmV0dXJuIGRldFxuICAgIH1cbiAgICByZXR1cm4gb3JpZW50YXRpb240RXhhY3QoYSxiLGMsZClcbiAgfVxuXVxuXG5mdW5jdGlvbiBzbG93T3JpZW50KGFyZ3MpIHtcbiAgdmFyIHByb2MgPSBDQUNIRURbYXJncy5sZW5ndGhdXG4gIGlmKCFwcm9jKSB7XG4gICAgcHJvYyA9IENBQ0hFRFthcmdzLmxlbmd0aF0gPSBvcmllbnRhdGlvbihhcmdzLmxlbmd0aClcbiAgfVxuICByZXR1cm4gcHJvYy5hcHBseSh1bmRlZmluZWQsIGFyZ3MpXG59XG5cbmZ1bmN0aW9uIGdlbmVyYXRlT3JpZW50YXRpb25Qcm9jKCkge1xuICB3aGlsZShDQUNIRUQubGVuZ3RoIDw9IE5VTV9FWFBBTkQpIHtcbiAgICBDQUNIRUQucHVzaChvcmllbnRhdGlvbihDQUNIRUQubGVuZ3RoKSlcbiAgfVxuICB2YXIgYXJncyA9IFtdXG4gIHZhciBwcm9jQXJncyA9IFtcInNsb3dcIl1cbiAgZm9yKHZhciBpPTA7IGk8PU5VTV9FWFBBTkQ7ICsraSkge1xuICAgIGFyZ3MucHVzaChcImFcIiArIGkpXG4gICAgcHJvY0FyZ3MucHVzaChcIm9cIiArIGkpXG4gIH1cbiAgdmFyIGNvZGUgPSBbXG4gICAgXCJmdW5jdGlvbiBnZXRPcmllbnRhdGlvbihcIiwgYXJncy5qb2luKCksIFwiKXtzd2l0Y2goYXJndW1lbnRzLmxlbmd0aCl7Y2FzZSAwOmNhc2UgMTpyZXR1cm4gMDtcIlxuICBdXG4gIGZvcih2YXIgaT0yOyBpPD1OVU1fRVhQQU5EOyArK2kpIHtcbiAgICBjb2RlLnB1c2goXCJjYXNlIFwiLCBpLCBcIjpyZXR1cm4gb1wiLCBpLCBcIihcIiwgYXJncy5zbGljZSgwLCBpKS5qb2luKCksIFwiKTtcIilcbiAgfVxuICBjb2RlLnB1c2goXCJ9dmFyIHM9bmV3IEFycmF5KGFyZ3VtZW50cy5sZW5ndGgpO2Zvcih2YXIgaT0wO2k8YXJndW1lbnRzLmxlbmd0aDsrK2kpe3NbaV09YXJndW1lbnRzW2ldfTtyZXR1cm4gc2xvdyhzKTt9cmV0dXJuIGdldE9yaWVudGF0aW9uXCIpXG4gIHByb2NBcmdzLnB1c2goY29kZS5qb2luKFwiXCIpKVxuXG4gIHZhciBwcm9jID0gRnVuY3Rpb24uYXBwbHkodW5kZWZpbmVkLCBwcm9jQXJncylcbiAgbW9kdWxlLmV4cG9ydHMgPSBwcm9jLmFwcGx5KHVuZGVmaW5lZCwgW3Nsb3dPcmllbnRdLmNvbmNhdChDQUNIRUQpKVxuICBmb3IodmFyIGk9MDsgaTw9TlVNX0VYUEFORDsgKytpKSB7XG4gICAgbW9kdWxlLmV4cG9ydHNbaV0gPSBDQUNIRURbaV1cbiAgfVxufVxuXG5nZW5lcmF0ZU9yaWVudGF0aW9uUHJvYygpIiwiLyoqXG4gKiBCaXQgdHdpZGRsaW5nIGhhY2tzIGZvciBKYXZhU2NyaXB0LlxuICpcbiAqIEF1dGhvcjogTWlrb2xhIEx5c2Vua29cbiAqXG4gKiBQb3J0ZWQgZnJvbSBTdGFuZm9yZCBiaXQgdHdpZGRsaW5nIGhhY2sgbGlicmFyeTpcbiAqICAgIGh0dHA6Ly9ncmFwaGljcy5zdGFuZm9yZC5lZHUvfnNlYW5kZXIvYml0aGFja3MuaHRtbFxuICovXG5cblwidXNlIHN0cmljdFwiOyBcInVzZSByZXN0cmljdFwiO1xuXG4vL051bWJlciBvZiBiaXRzIGluIGFuIGludGVnZXJcbnZhciBJTlRfQklUUyA9IDMyO1xuXG4vL0NvbnN0YW50c1xuZXhwb3J0cy5JTlRfQklUUyAgPSBJTlRfQklUUztcbmV4cG9ydHMuSU5UX01BWCAgID0gIDB4N2ZmZmZmZmY7XG5leHBvcnRzLklOVF9NSU4gICA9IC0xPDwoSU5UX0JJVFMtMSk7XG5cbi8vUmV0dXJucyAtMSwgMCwgKzEgZGVwZW5kaW5nIG9uIHNpZ24gb2YgeFxuZXhwb3J0cy5zaWduID0gZnVuY3Rpb24odikge1xuICByZXR1cm4gKHYgPiAwKSAtICh2IDwgMCk7XG59XG5cbi8vQ29tcHV0ZXMgYWJzb2x1dGUgdmFsdWUgb2YgaW50ZWdlclxuZXhwb3J0cy5hYnMgPSBmdW5jdGlvbih2KSB7XG4gIHZhciBtYXNrID0gdiA+PiAoSU5UX0JJVFMtMSk7XG4gIHJldHVybiAodiBeIG1hc2spIC0gbWFzaztcbn1cblxuLy9Db21wdXRlcyBtaW5pbXVtIG9mIGludGVnZXJzIHggYW5kIHlcbmV4cG9ydHMubWluID0gZnVuY3Rpb24oeCwgeSkge1xuICByZXR1cm4geSBeICgoeCBeIHkpICYgLSh4IDwgeSkpO1xufVxuXG4vL0NvbXB1dGVzIG1heGltdW0gb2YgaW50ZWdlcnMgeCBhbmQgeVxuZXhwb3J0cy5tYXggPSBmdW5jdGlvbih4LCB5KSB7XG4gIHJldHVybiB4IF4gKCh4IF4geSkgJiAtKHggPCB5KSk7XG59XG5cbi8vQ2hlY2tzIGlmIGEgbnVtYmVyIGlzIGEgcG93ZXIgb2YgdHdvXG5leHBvcnRzLmlzUG93MiA9IGZ1bmN0aW9uKHYpIHtcbiAgcmV0dXJuICEodiAmICh2LTEpKSAmJiAoISF2KTtcbn1cblxuLy9Db21wdXRlcyBsb2cgYmFzZSAyIG9mIHZcbmV4cG9ydHMubG9nMiA9IGZ1bmN0aW9uKHYpIHtcbiAgdmFyIHIsIHNoaWZ0O1xuICByID0gICAgICh2ID4gMHhGRkZGKSA8PCA0OyB2ID4+Pj0gcjtcbiAgc2hpZnQgPSAodiA+IDB4RkYgICkgPDwgMzsgdiA+Pj49IHNoaWZ0OyByIHw9IHNoaWZ0O1xuICBzaGlmdCA9ICh2ID4gMHhGICAgKSA8PCAyOyB2ID4+Pj0gc2hpZnQ7IHIgfD0gc2hpZnQ7XG4gIHNoaWZ0ID0gKHYgPiAweDMgICApIDw8IDE7IHYgPj4+PSBzaGlmdDsgciB8PSBzaGlmdDtcbiAgcmV0dXJuIHIgfCAodiA+PiAxKTtcbn1cblxuLy9Db21wdXRlcyBsb2cgYmFzZSAxMCBvZiB2XG5leHBvcnRzLmxvZzEwID0gZnVuY3Rpb24odikge1xuICByZXR1cm4gICh2ID49IDEwMDAwMDAwMDApID8gOSA6ICh2ID49IDEwMDAwMDAwMCkgPyA4IDogKHYgPj0gMTAwMDAwMDApID8gNyA6XG4gICAgICAgICAgKHYgPj0gMTAwMDAwMCkgPyA2IDogKHYgPj0gMTAwMDAwKSA/IDUgOiAodiA+PSAxMDAwMCkgPyA0IDpcbiAgICAgICAgICAodiA+PSAxMDAwKSA/IDMgOiAodiA+PSAxMDApID8gMiA6ICh2ID49IDEwKSA/IDEgOiAwO1xufVxuXG4vL0NvdW50cyBudW1iZXIgb2YgYml0c1xuZXhwb3J0cy5wb3BDb3VudCA9IGZ1bmN0aW9uKHYpIHtcbiAgdiA9IHYgLSAoKHYgPj4+IDEpICYgMHg1NTU1NTU1NSk7XG4gIHYgPSAodiAmIDB4MzMzMzMzMzMpICsgKCh2ID4+PiAyKSAmIDB4MzMzMzMzMzMpO1xuICByZXR1cm4gKCh2ICsgKHYgPj4+IDQpICYgMHhGMEYwRjBGKSAqIDB4MTAxMDEwMSkgPj4+IDI0O1xufVxuXG4vL0NvdW50cyBudW1iZXIgb2YgdHJhaWxpbmcgemVyb3NcbmZ1bmN0aW9uIGNvdW50VHJhaWxpbmdaZXJvcyh2KSB7XG4gIHZhciBjID0gMzI7XG4gIHYgJj0gLXY7XG4gIGlmICh2KSBjLS07XG4gIGlmICh2ICYgMHgwMDAwRkZGRikgYyAtPSAxNjtcbiAgaWYgKHYgJiAweDAwRkYwMEZGKSBjIC09IDg7XG4gIGlmICh2ICYgMHgwRjBGMEYwRikgYyAtPSA0O1xuICBpZiAodiAmIDB4MzMzMzMzMzMpIGMgLT0gMjtcbiAgaWYgKHYgJiAweDU1NTU1NTU1KSBjIC09IDE7XG4gIHJldHVybiBjO1xufVxuZXhwb3J0cy5jb3VudFRyYWlsaW5nWmVyb3MgPSBjb3VudFRyYWlsaW5nWmVyb3M7XG5cbi8vUm91bmRzIHRvIG5leHQgcG93ZXIgb2YgMlxuZXhwb3J0cy5uZXh0UG93MiA9IGZ1bmN0aW9uKHYpIHtcbiAgdiArPSB2ID09PSAwO1xuICAtLXY7XG4gIHYgfD0gdiA+Pj4gMTtcbiAgdiB8PSB2ID4+PiAyO1xuICB2IHw9IHYgPj4+IDQ7XG4gIHYgfD0gdiA+Pj4gODtcbiAgdiB8PSB2ID4+PiAxNjtcbiAgcmV0dXJuIHYgKyAxO1xufVxuXG4vL1JvdW5kcyBkb3duIHRvIHByZXZpb3VzIHBvd2VyIG9mIDJcbmV4cG9ydHMucHJldlBvdzIgPSBmdW5jdGlvbih2KSB7XG4gIHYgfD0gdiA+Pj4gMTtcbiAgdiB8PSB2ID4+PiAyO1xuICB2IHw9IHYgPj4+IDQ7XG4gIHYgfD0gdiA+Pj4gODtcbiAgdiB8PSB2ID4+PiAxNjtcbiAgcmV0dXJuIHYgLSAodj4+PjEpO1xufVxuXG4vL0NvbXB1dGVzIHBhcml0eSBvZiB3b3JkXG5leHBvcnRzLnBhcml0eSA9IGZ1bmN0aW9uKHYpIHtcbiAgdiBePSB2ID4+PiAxNjtcbiAgdiBePSB2ID4+PiA4O1xuICB2IF49IHYgPj4+IDQ7XG4gIHYgJj0gMHhmO1xuICByZXR1cm4gKDB4Njk5NiA+Pj4gdikgJiAxO1xufVxuXG52YXIgUkVWRVJTRV9UQUJMRSA9IG5ldyBBcnJheSgyNTYpO1xuXG4oZnVuY3Rpb24odGFiKSB7XG4gIGZvcih2YXIgaT0wOyBpPDI1NjsgKytpKSB7XG4gICAgdmFyIHYgPSBpLCByID0gaSwgcyA9IDc7XG4gICAgZm9yICh2ID4+Pj0gMTsgdjsgdiA+Pj49IDEpIHtcbiAgICAgIHIgPDw9IDE7XG4gICAgICByIHw9IHYgJiAxO1xuICAgICAgLS1zO1xuICAgIH1cbiAgICB0YWJbaV0gPSAociA8PCBzKSAmIDB4ZmY7XG4gIH1cbn0pKFJFVkVSU0VfVEFCTEUpO1xuXG4vL1JldmVyc2UgYml0cyBpbiBhIDMyIGJpdCB3b3JkXG5leHBvcnRzLnJldmVyc2UgPSBmdW5jdGlvbih2KSB7XG4gIHJldHVybiAgKFJFVkVSU0VfVEFCTEVbIHYgICAgICAgICAmIDB4ZmZdIDw8IDI0KSB8XG4gICAgICAgICAgKFJFVkVSU0VfVEFCTEVbKHYgPj4+IDgpICAmIDB4ZmZdIDw8IDE2KSB8XG4gICAgICAgICAgKFJFVkVSU0VfVEFCTEVbKHYgPj4+IDE2KSAmIDB4ZmZdIDw8IDgpICB8XG4gICAgICAgICAgIFJFVkVSU0VfVEFCTEVbKHYgPj4+IDI0KSAmIDB4ZmZdO1xufVxuXG4vL0ludGVybGVhdmUgYml0cyBvZiAyIGNvb3JkaW5hdGVzIHdpdGggMTYgYml0cy4gIFVzZWZ1bCBmb3IgZmFzdCBxdWFkdHJlZSBjb2Rlc1xuZXhwb3J0cy5pbnRlcmxlYXZlMiA9IGZ1bmN0aW9uKHgsIHkpIHtcbiAgeCAmPSAweEZGRkY7XG4gIHggPSAoeCB8ICh4IDw8IDgpKSAmIDB4MDBGRjAwRkY7XG4gIHggPSAoeCB8ICh4IDw8IDQpKSAmIDB4MEYwRjBGMEY7XG4gIHggPSAoeCB8ICh4IDw8IDIpKSAmIDB4MzMzMzMzMzM7XG4gIHggPSAoeCB8ICh4IDw8IDEpKSAmIDB4NTU1NTU1NTU7XG5cbiAgeSAmPSAweEZGRkY7XG4gIHkgPSAoeSB8ICh5IDw8IDgpKSAmIDB4MDBGRjAwRkY7XG4gIHkgPSAoeSB8ICh5IDw8IDQpKSAmIDB4MEYwRjBGMEY7XG4gIHkgPSAoeSB8ICh5IDw8IDIpKSAmIDB4MzMzMzMzMzM7XG4gIHkgPSAoeSB8ICh5IDw8IDEpKSAmIDB4NTU1NTU1NTU7XG5cbiAgcmV0dXJuIHggfCAoeSA8PCAxKTtcbn1cblxuLy9FeHRyYWN0cyB0aGUgbnRoIGludGVybGVhdmVkIGNvbXBvbmVudFxuZXhwb3J0cy5kZWludGVybGVhdmUyID0gZnVuY3Rpb24odiwgbikge1xuICB2ID0gKHYgPj4+IG4pICYgMHg1NTU1NTU1NTtcbiAgdiA9ICh2IHwgKHYgPj4+IDEpKSAgJiAweDMzMzMzMzMzO1xuICB2ID0gKHYgfCAodiA+Pj4gMikpICAmIDB4MEYwRjBGMEY7XG4gIHYgPSAodiB8ICh2ID4+PiA0KSkgICYgMHgwMEZGMDBGRjtcbiAgdiA9ICh2IHwgKHYgPj4+IDE2KSkgJiAweDAwMEZGRkY7XG4gIHJldHVybiAodiA8PCAxNikgPj4gMTY7XG59XG5cblxuLy9JbnRlcmxlYXZlIGJpdHMgb2YgMyBjb29yZGluYXRlcywgZWFjaCB3aXRoIDEwIGJpdHMuICBVc2VmdWwgZm9yIGZhc3Qgb2N0cmVlIGNvZGVzXG5leHBvcnRzLmludGVybGVhdmUzID0gZnVuY3Rpb24oeCwgeSwgeikge1xuICB4ICY9IDB4M0ZGO1xuICB4ICA9ICh4IHwgKHg8PDE2KSkgJiA0Mjc4MTkwMzM1O1xuICB4ICA9ICh4IHwgKHg8PDgpKSAgJiAyNTE3MTk2OTU7XG4gIHggID0gKHggfCAoeDw8NCkpICAmIDMyNzIzNTYwMzU7XG4gIHggID0gKHggfCAoeDw8MikpICAmIDEyMjcxMzM1MTM7XG5cbiAgeSAmPSAweDNGRjtcbiAgeSAgPSAoeSB8ICh5PDwxNikpICYgNDI3ODE5MDMzNTtcbiAgeSAgPSAoeSB8ICh5PDw4KSkgICYgMjUxNzE5Njk1O1xuICB5ICA9ICh5IHwgKHk8PDQpKSAgJiAzMjcyMzU2MDM1O1xuICB5ICA9ICh5IHwgKHk8PDIpKSAgJiAxMjI3MTMzNTEzO1xuICB4IHw9ICh5IDw8IDEpO1xuICBcbiAgeiAmPSAweDNGRjtcbiAgeiAgPSAoeiB8ICh6PDwxNikpICYgNDI3ODE5MDMzNTtcbiAgeiAgPSAoeiB8ICh6PDw4KSkgICYgMjUxNzE5Njk1O1xuICB6ICA9ICh6IHwgKHo8PDQpKSAgJiAzMjcyMzU2MDM1O1xuICB6ICA9ICh6IHwgKHo8PDIpKSAgJiAxMjI3MTMzNTEzO1xuICBcbiAgcmV0dXJuIHggfCAoeiA8PCAyKTtcbn1cblxuLy9FeHRyYWN0cyBudGggaW50ZXJsZWF2ZWQgY29tcG9uZW50IG9mIGEgMy10dXBsZVxuZXhwb3J0cy5kZWludGVybGVhdmUzID0gZnVuY3Rpb24odiwgbikge1xuICB2ID0gKHYgPj4+IG4pICAgICAgICYgMTIyNzEzMzUxMztcbiAgdiA9ICh2IHwgKHY+Pj4yKSkgICAmIDMyNzIzNTYwMzU7XG4gIHYgPSAodiB8ICh2Pj4+NCkpICAgJiAyNTE3MTk2OTU7XG4gIHYgPSAodiB8ICh2Pj4+OCkpICAgJiA0Mjc4MTkwMzM1O1xuICB2ID0gKHYgfCAodj4+PjE2KSkgICYgMHgzRkY7XG4gIHJldHVybiAodjw8MjIpPj4yMjtcbn1cblxuLy9Db21wdXRlcyBuZXh0IGNvbWJpbmF0aW9uIGluIGNvbGV4aWNvZ3JhcGhpYyBvcmRlciAodGhpcyBpcyBtaXN0YWtlbmx5IGNhbGxlZCBuZXh0UGVybXV0YXRpb24gb24gdGhlIGJpdCB0d2lkZGxpbmcgaGFja3MgcGFnZSlcbmV4cG9ydHMubmV4dENvbWJpbmF0aW9uID0gZnVuY3Rpb24odikge1xuICB2YXIgdCA9IHYgfCAodiAtIDEpO1xuICByZXR1cm4gKHQgKyAxKSB8ICgoKH50ICYgLX50KSAtIDEpID4+PiAoY291bnRUcmFpbGluZ1plcm9zKHYpICsgMSkpO1xufVxuXG4iLCJcInVzZSBzdHJpY3RcIjsgXCJ1c2UgcmVzdHJpY3RcIjtcblxubW9kdWxlLmV4cG9ydHMgPSBVbmlvbkZpbmQ7XG5cbmZ1bmN0aW9uIFVuaW9uRmluZChjb3VudCkge1xuICB0aGlzLnJvb3RzID0gbmV3IEFycmF5KGNvdW50KTtcbiAgdGhpcy5yYW5rcyA9IG5ldyBBcnJheShjb3VudCk7XG4gIFxuICBmb3IodmFyIGk9MDsgaTxjb3VudDsgKytpKSB7XG4gICAgdGhpcy5yb290c1tpXSA9IGk7XG4gICAgdGhpcy5yYW5rc1tpXSA9IDA7XG4gIH1cbn1cblxudmFyIHByb3RvID0gVW5pb25GaW5kLnByb3RvdHlwZVxuXG5PYmplY3QuZGVmaW5lUHJvcGVydHkocHJvdG8sIFwibGVuZ3RoXCIsIHtcbiAgXCJnZXRcIjogZnVuY3Rpb24oKSB7XG4gICAgcmV0dXJuIHRoaXMucm9vdHMubGVuZ3RoXG4gIH1cbn0pXG5cbnByb3RvLm1ha2VTZXQgPSBmdW5jdGlvbigpIHtcbiAgdmFyIG4gPSB0aGlzLnJvb3RzLmxlbmd0aDtcbiAgdGhpcy5yb290cy5wdXNoKG4pO1xuICB0aGlzLnJhbmtzLnB1c2goMCk7XG4gIHJldHVybiBuO1xufVxuXG5wcm90by5maW5kID0gZnVuY3Rpb24oeCkge1xuICB2YXIgcm9vdHMgPSB0aGlzLnJvb3RzO1xuICB3aGlsZShyb290c1t4XSAhPT0geCkge1xuICAgIHZhciB5ID0gcm9vdHNbeF07XG4gICAgcm9vdHNbeF0gPSByb290c1t5XTtcbiAgICB4ID0geTtcbiAgfVxuICByZXR1cm4geDtcbn1cblxucHJvdG8ubGluayA9IGZ1bmN0aW9uKHgsIHkpIHtcbiAgdmFyIHhyID0gdGhpcy5maW5kKHgpXG4gICAgLCB5ciA9IHRoaXMuZmluZCh5KTtcbiAgaWYoeHIgPT09IHlyKSB7XG4gICAgcmV0dXJuO1xuICB9XG4gIHZhciByYW5rcyA9IHRoaXMucmFua3NcbiAgICAsIHJvb3RzID0gdGhpcy5yb290c1xuICAgICwgeGQgICAgPSByYW5rc1t4cl1cbiAgICAsIHlkICAgID0gcmFua3NbeXJdO1xuICBpZih4ZCA8IHlkKSB7XG4gICAgcm9vdHNbeHJdID0geXI7XG4gIH0gZWxzZSBpZih5ZCA8IHhkKSB7XG4gICAgcm9vdHNbeXJdID0geHI7XG4gIH0gZWxzZSB7XG4gICAgcm9vdHNbeXJdID0geHI7XG4gICAgKytyYW5rc1t4cl07XG4gIH1cbn0iLCJcInVzZSBzdHJpY3RcIjsgXCJ1c2UgcmVzdHJpY3RcIjtcblxudmFyIGJpdHMgICAgICA9IHJlcXVpcmUoXCJiaXQtdHdpZGRsZVwiKVxuICAsIFVuaW9uRmluZCA9IHJlcXVpcmUoXCJ1bmlvbi1maW5kXCIpXG5cbi8vUmV0dXJucyB0aGUgZGltZW5zaW9uIG9mIGEgY2VsbCBjb21wbGV4XG5mdW5jdGlvbiBkaW1lbnNpb24oY2VsbHMpIHtcbiAgdmFyIGQgPSAwXG4gICAgLCBtYXggPSBNYXRoLm1heFxuICBmb3IodmFyIGk9MCwgaWw9Y2VsbHMubGVuZ3RoOyBpPGlsOyArK2kpIHtcbiAgICBkID0gbWF4KGQsIGNlbGxzW2ldLmxlbmd0aClcbiAgfVxuICByZXR1cm4gZC0xXG59XG5leHBvcnRzLmRpbWVuc2lvbiA9IGRpbWVuc2lvblxuXG4vL0NvdW50cyB0aGUgbnVtYmVyIG9mIHZlcnRpY2VzIGluIGZhY2VzXG5mdW5jdGlvbiBjb3VudFZlcnRpY2VzKGNlbGxzKSB7XG4gIHZhciB2YyA9IC0xXG4gICAgLCBtYXggPSBNYXRoLm1heFxuICBmb3IodmFyIGk9MCwgaWw9Y2VsbHMubGVuZ3RoOyBpPGlsOyArK2kpIHtcbiAgICB2YXIgYyA9IGNlbGxzW2ldXG4gICAgZm9yKHZhciBqPTAsIGpsPWMubGVuZ3RoOyBqPGpsOyArK2opIHtcbiAgICAgIHZjID0gbWF4KHZjLCBjW2pdKVxuICAgIH1cbiAgfVxuICByZXR1cm4gdmMrMVxufVxuZXhwb3J0cy5jb3VudFZlcnRpY2VzID0gY291bnRWZXJ0aWNlc1xuXG4vL1JldHVybnMgYSBkZWVwIGNvcHkgb2YgY2VsbHNcbmZ1bmN0aW9uIGNsb25lQ2VsbHMoY2VsbHMpIHtcbiAgdmFyIG5jZWxscyA9IG5ldyBBcnJheShjZWxscy5sZW5ndGgpXG4gIGZvcih2YXIgaT0wLCBpbD1jZWxscy5sZW5ndGg7IGk8aWw7ICsraSkge1xuICAgIG5jZWxsc1tpXSA9IGNlbGxzW2ldLnNsaWNlKDApXG4gIH1cbiAgcmV0dXJuIG5jZWxsc1xufVxuZXhwb3J0cy5jbG9uZUNlbGxzID0gY2xvbmVDZWxsc1xuXG4vL1JhbmtzIGEgcGFpciBvZiBjZWxscyB1cCB0byBwZXJtdXRhdGlvblxuZnVuY3Rpb24gY29tcGFyZUNlbGxzKGEsIGIpIHtcbiAgdmFyIG4gPSBhLmxlbmd0aFxuICAgICwgdCA9IGEubGVuZ3RoIC0gYi5sZW5ndGhcbiAgICAsIG1pbiA9IE1hdGgubWluXG4gIGlmKHQpIHtcbiAgICByZXR1cm4gdFxuICB9XG4gIHN3aXRjaChuKSB7XG4gICAgY2FzZSAwOlxuICAgICAgcmV0dXJuIDA7XG4gICAgY2FzZSAxOlxuICAgICAgcmV0dXJuIGFbMF0gLSBiWzBdO1xuICAgIGNhc2UgMjpcbiAgICAgIHZhciBkID0gYVswXSthWzFdLWJbMF0tYlsxXVxuICAgICAgaWYoZCkge1xuICAgICAgICByZXR1cm4gZFxuICAgICAgfVxuICAgICAgcmV0dXJuIG1pbihhWzBdLGFbMV0pIC0gbWluKGJbMF0sYlsxXSlcbiAgICBjYXNlIDM6XG4gICAgICB2YXIgbDEgPSBhWzBdK2FbMV1cbiAgICAgICAgLCBtMSA9IGJbMF0rYlsxXVxuICAgICAgZCA9IGwxK2FbMl0gLSAobTErYlsyXSlcbiAgICAgIGlmKGQpIHtcbiAgICAgICAgcmV0dXJuIGRcbiAgICAgIH1cbiAgICAgIHZhciBsMCA9IG1pbihhWzBdLCBhWzFdKVxuICAgICAgICAsIG0wID0gbWluKGJbMF0sIGJbMV0pXG4gICAgICAgICwgZCAgPSBtaW4obDAsIGFbMl0pIC0gbWluKG0wLCBiWzJdKVxuICAgICAgaWYoZCkge1xuICAgICAgICByZXR1cm4gZFxuICAgICAgfVxuICAgICAgcmV0dXJuIG1pbihsMCthWzJdLCBsMSkgLSBtaW4obTArYlsyXSwgbTEpXG4gICAgXG4gICAgLy9UT0RPOiBNYXliZSBvcHRpbWl6ZSBuPTQgYXMgd2VsbD9cbiAgICBcbiAgICBkZWZhdWx0OlxuICAgICAgdmFyIGFzID0gYS5zbGljZSgwKVxuICAgICAgYXMuc29ydCgpXG4gICAgICB2YXIgYnMgPSBiLnNsaWNlKDApXG4gICAgICBicy5zb3J0KClcbiAgICAgIGZvcih2YXIgaT0wOyBpPG47ICsraSkge1xuICAgICAgICB0ID0gYXNbaV0gLSBic1tpXVxuICAgICAgICBpZih0KSB7XG4gICAgICAgICAgcmV0dXJuIHRcbiAgICAgICAgfVxuICAgICAgfVxuICAgICAgcmV0dXJuIDBcbiAgfVxufVxuZXhwb3J0cy5jb21wYXJlQ2VsbHMgPSBjb21wYXJlQ2VsbHNcblxuZnVuY3Rpb24gY29tcGFyZVppcHBlZChhLCBiKSB7XG4gIHJldHVybiBjb21wYXJlQ2VsbHMoYVswXSwgYlswXSlcbn1cblxuLy9QdXRzIGEgY2VsbCBjb21wbGV4IGludG8gbm9ybWFsIG9yZGVyIGZvciB0aGUgcHVycG9zZXMgb2YgZmluZENlbGwgcXVlcmllc1xuZnVuY3Rpb24gbm9ybWFsaXplKGNlbGxzLCBhdHRyKSB7XG4gIGlmKGF0dHIpIHtcbiAgICB2YXIgbGVuID0gY2VsbHMubGVuZ3RoXG4gICAgdmFyIHppcHBlZCA9IG5ldyBBcnJheShsZW4pXG4gICAgZm9yKHZhciBpPTA7IGk8bGVuOyArK2kpIHtcbiAgICAgIHppcHBlZFtpXSA9IFtjZWxsc1tpXSwgYXR0cltpXV1cbiAgICB9XG4gICAgemlwcGVkLnNvcnQoY29tcGFyZVppcHBlZClcbiAgICBmb3IodmFyIGk9MDsgaTxsZW47ICsraSkge1xuICAgICAgY2VsbHNbaV0gPSB6aXBwZWRbaV1bMF1cbiAgICAgIGF0dHJbaV0gPSB6aXBwZWRbaV1bMV1cbiAgICB9XG4gICAgcmV0dXJuIGNlbGxzXG4gIH0gZWxzZSB7XG4gICAgY2VsbHMuc29ydChjb21wYXJlQ2VsbHMpXG4gICAgcmV0dXJuIGNlbGxzXG4gIH1cbn1cbmV4cG9ydHMubm9ybWFsaXplID0gbm9ybWFsaXplXG5cbi8vUmVtb3ZlcyBhbGwgZHVwbGljYXRlIGNlbGxzIGluIHRoZSBjb21wbGV4XG5mdW5jdGlvbiB1bmlxdWUoY2VsbHMpIHtcbiAgaWYoY2VsbHMubGVuZ3RoID09PSAwKSB7XG4gICAgcmV0dXJuIFtdXG4gIH1cbiAgdmFyIHB0ciA9IDFcbiAgICAsIGxlbiA9IGNlbGxzLmxlbmd0aFxuICBmb3IodmFyIGk9MTsgaTxsZW47ICsraSkge1xuICAgIHZhciBhID0gY2VsbHNbaV1cbiAgICBpZihjb21wYXJlQ2VsbHMoYSwgY2VsbHNbaS0xXSkpIHtcbiAgICAgIGlmKGkgPT09IHB0cikge1xuICAgICAgICBwdHIrK1xuICAgICAgICBjb250aW51ZVxuICAgICAgfVxuICAgICAgY2VsbHNbcHRyKytdID0gYVxuICAgIH1cbiAgfVxuICBjZWxscy5sZW5ndGggPSBwdHJcbiAgcmV0dXJuIGNlbGxzXG59XG5leHBvcnRzLnVuaXF1ZSA9IHVuaXF1ZTtcblxuLy9GaW5kcyBhIGNlbGwgaW4gYSBub3JtYWxpemVkIGNlbGwgY29tcGxleFxuZnVuY3Rpb24gZmluZENlbGwoY2VsbHMsIGMpIHtcbiAgdmFyIGxvID0gMFxuICAgICwgaGkgPSBjZWxscy5sZW5ndGgtMVxuICAgICwgciAgPSAtMVxuICB3aGlsZSAobG8gPD0gaGkpIHtcbiAgICB2YXIgbWlkID0gKGxvICsgaGkpID4+IDFcbiAgICAgICwgcyAgID0gY29tcGFyZUNlbGxzKGNlbGxzW21pZF0sIGMpXG4gICAgaWYocyA8PSAwKSB7XG4gICAgICBpZihzID09PSAwKSB7XG4gICAgICAgIHIgPSBtaWRcbiAgICAgIH1cbiAgICAgIGxvID0gbWlkICsgMVxuICAgIH0gZWxzZSBpZihzID4gMCkge1xuICAgICAgaGkgPSBtaWQgLSAxXG4gICAgfVxuICB9XG4gIHJldHVybiByXG59XG5leHBvcnRzLmZpbmRDZWxsID0gZmluZENlbGw7XG5cbi8vQnVpbGRzIGFuIGluZGV4IGZvciBhbiBuLWNlbGwuICBUaGlzIGlzIG1vcmUgZ2VuZXJhbCB0aGFuIGR1YWwsIGJ1dCBsZXNzIGVmZmljaWVudFxuZnVuY3Rpb24gaW5jaWRlbmNlKGZyb21fY2VsbHMsIHRvX2NlbGxzKSB7XG4gIHZhciBpbmRleCA9IG5ldyBBcnJheShmcm9tX2NlbGxzLmxlbmd0aClcbiAgZm9yKHZhciBpPTAsIGlsPWluZGV4Lmxlbmd0aDsgaTxpbDsgKytpKSB7XG4gICAgaW5kZXhbaV0gPSBbXVxuICB9XG4gIHZhciBiID0gW11cbiAgZm9yKHZhciBpPTAsIG49dG9fY2VsbHMubGVuZ3RoOyBpPG47ICsraSkge1xuICAgIHZhciBjID0gdG9fY2VsbHNbaV1cbiAgICB2YXIgY2wgPSBjLmxlbmd0aFxuICAgIGZvcih2YXIgaz0xLCBrbj0oMTw8Y2wpOyBrPGtuOyArK2spIHtcbiAgICAgIGIubGVuZ3RoID0gYml0cy5wb3BDb3VudChrKVxuICAgICAgdmFyIGwgPSAwXG4gICAgICBmb3IodmFyIGo9MDsgajxjbDsgKytqKSB7XG4gICAgICAgIGlmKGsgJiAoMTw8aikpIHtcbiAgICAgICAgICBiW2wrK10gPSBjW2pdXG4gICAgICAgIH1cbiAgICAgIH1cbiAgICAgIHZhciBpZHg9ZmluZENlbGwoZnJvbV9jZWxscywgYilcbiAgICAgIGlmKGlkeCA8IDApIHtcbiAgICAgICAgY29udGludWVcbiAgICAgIH1cbiAgICAgIHdoaWxlKHRydWUpIHtcbiAgICAgICAgaW5kZXhbaWR4KytdLnB1c2goaSlcbiAgICAgICAgaWYoaWR4ID49IGZyb21fY2VsbHMubGVuZ3RoIHx8IGNvbXBhcmVDZWxscyhmcm9tX2NlbGxzW2lkeF0sIGIpICE9PSAwKSB7XG4gICAgICAgICAgYnJlYWtcbiAgICAgICAgfVxuICAgICAgfVxuICAgIH1cbiAgfVxuICByZXR1cm4gaW5kZXhcbn1cbmV4cG9ydHMuaW5jaWRlbmNlID0gaW5jaWRlbmNlXG5cbi8vQ29tcHV0ZXMgdGhlIGR1YWwgb2YgdGhlIG1lc2guICBUaGlzIGlzIGJhc2ljYWxseSBhbiBvcHRpbWl6ZWQgdmVyc2lvbiBvZiBidWlsZEluZGV4IGZvciB0aGUgc2l0dWF0aW9uIHdoZXJlIGZyb21fY2VsbHMgaXMganVzdCB0aGUgbGlzdCBvZiB2ZXJ0aWNlc1xuZnVuY3Rpb24gZHVhbChjZWxscywgdmVydGV4X2NvdW50KSB7XG4gIGlmKCF2ZXJ0ZXhfY291bnQpIHtcbiAgICByZXR1cm4gaW5jaWRlbmNlKHVuaXF1ZShza2VsZXRvbihjZWxscywgMCkpLCBjZWxscywgMClcbiAgfVxuICB2YXIgcmVzID0gbmV3IEFycmF5KHZlcnRleF9jb3VudClcbiAgZm9yKHZhciBpPTA7IGk8dmVydGV4X2NvdW50OyArK2kpIHtcbiAgICByZXNbaV0gPSBbXVxuICB9XG4gIGZvcih2YXIgaT0wLCBsZW49Y2VsbHMubGVuZ3RoOyBpPGxlbjsgKytpKSB7XG4gICAgdmFyIGMgPSBjZWxsc1tpXVxuICAgIGZvcih2YXIgaj0wLCBjbD1jLmxlbmd0aDsgajxjbDsgKytqKSB7XG4gICAgICByZXNbY1tqXV0ucHVzaChpKVxuICAgIH1cbiAgfVxuICByZXR1cm4gcmVzXG59XG5leHBvcnRzLmR1YWwgPSBkdWFsXG5cbi8vRW51bWVyYXRlcyBhbGwgY2VsbHMgaW4gdGhlIGNvbXBsZXhcbmZ1bmN0aW9uIGV4cGxvZGUoY2VsbHMpIHtcbiAgdmFyIHJlc3VsdCA9IFtdXG4gIGZvcih2YXIgaT0wLCBpbD1jZWxscy5sZW5ndGg7IGk8aWw7ICsraSkge1xuICAgIHZhciBjID0gY2VsbHNbaV1cbiAgICAgICwgY2wgPSBjLmxlbmd0aHwwXG4gICAgZm9yKHZhciBqPTEsIGpsPSgxPDxjbCk7IGo8amw7ICsraikge1xuICAgICAgdmFyIGIgPSBbXVxuICAgICAgZm9yKHZhciBrPTA7IGs8Y2w7ICsraykge1xuICAgICAgICBpZigoaiA+Pj4gaykgJiAxKSB7XG4gICAgICAgICAgYi5wdXNoKGNba10pXG4gICAgICAgIH1cbiAgICAgIH1cbiAgICAgIHJlc3VsdC5wdXNoKGIpXG4gICAgfVxuICB9XG4gIHJldHVybiBub3JtYWxpemUocmVzdWx0KVxufVxuZXhwb3J0cy5leHBsb2RlID0gZXhwbG9kZVxuXG4vL0VudW1lcmF0ZXMgYWxsIG9mIHRoZSBuLWNlbGxzIG9mIGEgY2VsbCBjb21wbGV4XG5mdW5jdGlvbiBza2VsZXRvbihjZWxscywgbikge1xuICBpZihuIDwgMCkge1xuICAgIHJldHVybiBbXVxuICB9XG4gIHZhciByZXN1bHQgPSBbXVxuICAgICwgazAgICAgID0gKDE8PChuKzEpKS0xXG4gIGZvcih2YXIgaT0wOyBpPGNlbGxzLmxlbmd0aDsgKytpKSB7XG4gICAgdmFyIGMgPSBjZWxsc1tpXVxuICAgIGZvcih2YXIgaz1rMDsgazwoMTw8Yy5sZW5ndGgpOyBrPWJpdHMubmV4dENvbWJpbmF0aW9uKGspKSB7XG4gICAgICB2YXIgYiA9IG5ldyBBcnJheShuKzEpXG4gICAgICAgICwgbCA9IDBcbiAgICAgIGZvcih2YXIgaj0wOyBqPGMubGVuZ3RoOyArK2opIHtcbiAgICAgICAgaWYoayAmICgxPDxqKSkge1xuICAgICAgICAgIGJbbCsrXSA9IGNbal1cbiAgICAgICAgfVxuICAgICAgfVxuICAgICAgcmVzdWx0LnB1c2goYilcbiAgICB9XG4gIH1cbiAgcmV0dXJuIG5vcm1hbGl6ZShyZXN1bHQpXG59XG5leHBvcnRzLnNrZWxldG9uID0gc2tlbGV0b247XG5cbi8vQ29tcHV0ZXMgdGhlIGJvdW5kYXJ5IG9mIGFsbCBjZWxscywgZG9lcyBub3QgcmVtb3ZlIGR1cGxpY2F0ZXNcbmZ1bmN0aW9uIGJvdW5kYXJ5KGNlbGxzKSB7XG4gIHZhciByZXMgPSBbXVxuICBmb3IodmFyIGk9MCxpbD1jZWxscy5sZW5ndGg7IGk8aWw7ICsraSkge1xuICAgIHZhciBjID0gY2VsbHNbaV1cbiAgICBmb3IodmFyIGo9MCxjbD1jLmxlbmd0aDsgajxjbDsgKytqKSB7XG4gICAgICB2YXIgYiA9IG5ldyBBcnJheShjLmxlbmd0aC0xKVxuICAgICAgZm9yKHZhciBrPTAsIGw9MDsgazxjbDsgKytrKSB7XG4gICAgICAgIGlmKGsgIT09IGopIHtcbiAgICAgICAgICBiW2wrK10gPSBjW2tdXG4gICAgICAgIH1cbiAgICAgIH1cbiAgICAgIHJlcy5wdXNoKGIpXG4gICAgfVxuICB9XG4gIHJldHVybiBub3JtYWxpemUocmVzKVxufVxuZXhwb3J0cy5ib3VuZGFyeSA9IGJvdW5kYXJ5O1xuXG4vL0NvbXB1dGVzIGNvbm5lY3RlZCBjb21wb25lbnRzIGZvciBhIGRlbnNlIGNlbGwgY29tcGxleFxuZnVuY3Rpb24gY29ubmVjdGVkQ29tcG9uZW50c19kZW5zZShjZWxscywgdmVydGV4X2NvdW50KSB7XG4gIHZhciBsYWJlbHMgPSBuZXcgVW5pb25GaW5kKHZlcnRleF9jb3VudClcbiAgZm9yKHZhciBpPTA7IGk8Y2VsbHMubGVuZ3RoOyArK2kpIHtcbiAgICB2YXIgYyA9IGNlbGxzW2ldXG4gICAgZm9yKHZhciBqPTA7IGo8Yy5sZW5ndGg7ICsraikge1xuICAgICAgZm9yKHZhciBrPWorMTsgazxjLmxlbmd0aDsgKytrKSB7XG4gICAgICAgIGxhYmVscy5saW5rKGNbal0sIGNba10pXG4gICAgICB9XG4gICAgfVxuICB9XG4gIHZhciBjb21wb25lbnRzID0gW11cbiAgICAsIGNvbXBvbmVudF9sYWJlbHMgPSBsYWJlbHMucmFua3NcbiAgZm9yKHZhciBpPTA7IGk8Y29tcG9uZW50X2xhYmVscy5sZW5ndGg7ICsraSkge1xuICAgIGNvbXBvbmVudF9sYWJlbHNbaV0gPSAtMVxuICB9XG4gIGZvcih2YXIgaT0wOyBpPGNlbGxzLmxlbmd0aDsgKytpKSB7XG4gICAgdmFyIGwgPSBsYWJlbHMuZmluZChjZWxsc1tpXVswXSlcbiAgICBpZihjb21wb25lbnRfbGFiZWxzW2xdIDwgMCkge1xuICAgICAgY29tcG9uZW50X2xhYmVsc1tsXSA9IGNvbXBvbmVudHMubGVuZ3RoXG4gICAgICBjb21wb25lbnRzLnB1c2goW2NlbGxzW2ldLnNsaWNlKDApXSlcbiAgICB9IGVsc2Uge1xuICAgICAgY29tcG9uZW50c1tjb21wb25lbnRfbGFiZWxzW2xdXS5wdXNoKGNlbGxzW2ldLnNsaWNlKDApKVxuICAgIH1cbiAgfVxuICByZXR1cm4gY29tcG9uZW50c1xufVxuXG4vL0NvbXB1dGVzIGNvbm5lY3RlZCBjb21wb25lbnRzIGZvciBhIHNwYXJzZSBncmFwaFxuZnVuY3Rpb24gY29ubmVjdGVkQ29tcG9uZW50c19zcGFyc2UoY2VsbHMpIHtcbiAgdmFyIHZlcnRpY2VzICA9IHVuaXF1ZShub3JtYWxpemUoc2tlbGV0b24oY2VsbHMsIDApKSlcbiAgICAsIGxhYmVscyAgICA9IG5ldyBVbmlvbkZpbmQodmVydGljZXMubGVuZ3RoKVxuICBmb3IodmFyIGk9MDsgaTxjZWxscy5sZW5ndGg7ICsraSkge1xuICAgIHZhciBjID0gY2VsbHNbaV1cbiAgICBmb3IodmFyIGo9MDsgajxjLmxlbmd0aDsgKytqKSB7XG4gICAgICB2YXIgdmogPSBmaW5kQ2VsbCh2ZXJ0aWNlcywgW2Nbal1dKVxuICAgICAgZm9yKHZhciBrPWorMTsgazxjLmxlbmd0aDsgKytrKSB7XG4gICAgICAgIGxhYmVscy5saW5rKHZqLCBmaW5kQ2VsbCh2ZXJ0aWNlcywgW2Nba11dKSlcbiAgICAgIH1cbiAgICB9XG4gIH1cbiAgdmFyIGNvbXBvbmVudHMgICAgICAgID0gW11cbiAgICAsIGNvbXBvbmVudF9sYWJlbHMgID0gbGFiZWxzLnJhbmtzXG4gIGZvcih2YXIgaT0wOyBpPGNvbXBvbmVudF9sYWJlbHMubGVuZ3RoOyArK2kpIHtcbiAgICBjb21wb25lbnRfbGFiZWxzW2ldID0gLTFcbiAgfVxuICBmb3IodmFyIGk9MDsgaTxjZWxscy5sZW5ndGg7ICsraSkge1xuICAgIHZhciBsID0gbGFiZWxzLmZpbmQoZmluZENlbGwodmVydGljZXMsIFtjZWxsc1tpXVswXV0pKTtcbiAgICBpZihjb21wb25lbnRfbGFiZWxzW2xdIDwgMCkge1xuICAgICAgY29tcG9uZW50X2xhYmVsc1tsXSA9IGNvbXBvbmVudHMubGVuZ3RoXG4gICAgICBjb21wb25lbnRzLnB1c2goW2NlbGxzW2ldLnNsaWNlKDApXSlcbiAgICB9IGVsc2Uge1xuICAgICAgY29tcG9uZW50c1tjb21wb25lbnRfbGFiZWxzW2xdXS5wdXNoKGNlbGxzW2ldLnNsaWNlKDApKVxuICAgIH1cbiAgfVxuICByZXR1cm4gY29tcG9uZW50c1xufVxuXG4vL0NvbXB1dGVzIGNvbm5lY3RlZCBjb21wb25lbnRzIGZvciBhIGNlbGwgY29tcGxleFxuZnVuY3Rpb24gY29ubmVjdGVkQ29tcG9uZW50cyhjZWxscywgdmVydGV4X2NvdW50KSB7XG4gIGlmKHZlcnRleF9jb3VudCkge1xuICAgIHJldHVybiBjb25uZWN0ZWRDb21wb25lbnRzX2RlbnNlKGNlbGxzLCB2ZXJ0ZXhfY291bnQpXG4gIH1cbiAgcmV0dXJuIGNvbm5lY3RlZENvbXBvbmVudHNfc3BhcnNlKGNlbGxzKVxufVxuZXhwb3J0cy5jb25uZWN0ZWRDb21wb25lbnRzID0gY29ubmVjdGVkQ29tcG9uZW50c1xuIiwiXCJ1c2Ugc3RyaWN0XCJcblxuZnVuY3Rpb24gdW5pcXVlX3ByZWQobGlzdCwgY29tcGFyZSkge1xuICB2YXIgcHRyID0gMVxuICAgICwgbGVuID0gbGlzdC5sZW5ndGhcbiAgICAsIGE9bGlzdFswXSwgYj1saXN0WzBdXG4gIGZvcih2YXIgaT0xOyBpPGxlbjsgKytpKSB7XG4gICAgYiA9IGFcbiAgICBhID0gbGlzdFtpXVxuICAgIGlmKGNvbXBhcmUoYSwgYikpIHtcbiAgICAgIGlmKGkgPT09IHB0cikge1xuICAgICAgICBwdHIrK1xuICAgICAgICBjb250aW51ZVxuICAgICAgfVxuICAgICAgbGlzdFtwdHIrK10gPSBhXG4gICAgfVxuICB9XG4gIGxpc3QubGVuZ3RoID0gcHRyXG4gIHJldHVybiBsaXN0XG59XG5cbmZ1bmN0aW9uIHVuaXF1ZV9lcShsaXN0KSB7XG4gIHZhciBwdHIgPSAxXG4gICAgLCBsZW4gPSBsaXN0Lmxlbmd0aFxuICAgICwgYT1saXN0WzBdLCBiID0gbGlzdFswXVxuICBmb3IodmFyIGk9MTsgaTxsZW47ICsraSwgYj1hKSB7XG4gICAgYiA9IGFcbiAgICBhID0gbGlzdFtpXVxuICAgIGlmKGEgIT09IGIpIHtcbiAgICAgIGlmKGkgPT09IHB0cikge1xuICAgICAgICBwdHIrK1xuICAgICAgICBjb250aW51ZVxuICAgICAgfVxuICAgICAgbGlzdFtwdHIrK10gPSBhXG4gICAgfVxuICB9XG4gIGxpc3QubGVuZ3RoID0gcHRyXG4gIHJldHVybiBsaXN0XG59XG5cbmZ1bmN0aW9uIHVuaXF1ZShsaXN0LCBjb21wYXJlLCBzb3J0ZWQpIHtcbiAgaWYobGlzdC5sZW5ndGggPT09IDApIHtcbiAgICByZXR1cm4gbGlzdFxuICB9XG4gIGlmKGNvbXBhcmUpIHtcbiAgICBpZighc29ydGVkKSB7XG4gICAgICBsaXN0LnNvcnQoY29tcGFyZSlcbiAgICB9XG4gICAgcmV0dXJuIHVuaXF1ZV9wcmVkKGxpc3QsIGNvbXBhcmUpXG4gIH1cbiAgaWYoIXNvcnRlZCkge1xuICAgIGxpc3Quc29ydCgpXG4gIH1cbiAgcmV0dXJuIHVuaXF1ZV9lcShsaXN0KVxufVxuXG5tb2R1bGUuZXhwb3J0cyA9IHVuaXF1ZVxuIiwiXCJ1c2Ugc3RyaWN0XCJcblxudmFyIGNoID0gcmVxdWlyZShcImluY3JlbWVudGFsLWNvbnZleC1odWxsXCIpXG52YXIgdW5pcSA9IHJlcXVpcmUoXCJ1bmlxXCIpXG5cbm1vZHVsZS5leHBvcnRzID0gdHJpYW5ndWxhdGVcblxuZnVuY3Rpb24gTGlmdGVkUG9pbnQocCwgaSkge1xuICB0aGlzLnBvaW50ID0gcFxuICB0aGlzLmluZGV4ID0gaVxufVxuXG5mdW5jdGlvbiBjb21wYXJlTGlmdGVkKGEsIGIpIHtcbiAgdmFyIGFwID0gYS5wb2ludFxuICB2YXIgYnAgPSBiLnBvaW50XG4gIHZhciBkID0gYXAubGVuZ3RoXG4gIGZvcih2YXIgaT0wOyBpPGQ7ICsraSkge1xuICAgIHZhciBzID0gYnBbaV0gLSBhcFtpXVxuICAgIGlmKHMpIHtcbiAgICAgIHJldHVybiBzXG4gICAgfVxuICB9XG4gIHJldHVybiAwXG59XG5cbmZ1bmN0aW9uIHRyaWFuZ3VsYXRlMUQobiwgcG9pbnRzLCBpbmNsdWRlUG9pbnRBdEluZmluaXR5KSB7XG4gIGlmKG4gPT09IDEpIHtcbiAgICBpZihpbmNsdWRlUG9pbnRBdEluZmluaXR5KSB7XG4gICAgICByZXR1cm4gWyBbLTEsIDBdIF1cbiAgICB9IGVsc2Uge1xuICAgICAgcmV0dXJuIFtdXG4gICAgfVxuICB9XG4gIHZhciBsaWZ0ZWQgPSBwb2ludHMubWFwKGZ1bmN0aW9uKHAsIGkpIHtcbiAgICByZXR1cm4gWyBwWzBdLCBpIF1cbiAgfSlcbiAgbGlmdGVkLnNvcnQoZnVuY3Rpb24oYSxiKSB7XG4gICAgcmV0dXJuIGFbMF0gLSBiWzBdXG4gIH0pXG4gIHZhciBjZWxscyA9IG5ldyBBcnJheShuIC0gMSlcbiAgZm9yKHZhciBpPTE7IGk8bjsgKytpKSB7XG4gICAgdmFyIGEgPSBsaWZ0ZWRbaS0xXVxuICAgIHZhciBiID0gbGlmdGVkW2ldXG4gICAgY2VsbHNbaS0xXSA9IFsgYVsxXSwgYlsxXSBdXG4gIH1cbiAgaWYoaW5jbHVkZVBvaW50QXRJbmZpbml0eSkge1xuICAgIGNlbGxzLnB1c2goXG4gICAgICBbIC0xLCBjZWxsc1swXVsxXSwgXSxcbiAgICAgIFsgY2VsbHNbbi0xXVsxXSwgLTEgXSlcbiAgfVxuICByZXR1cm4gY2VsbHNcbn1cblxuZnVuY3Rpb24gdHJpYW5ndWxhdGUocG9pbnRzLCBpbmNsdWRlUG9pbnRBdEluZmluaXR5KSB7XG4gIHZhciBuID0gcG9pbnRzLmxlbmd0aFxuICBpZihuID09PSAwKSB7XG4gICAgcmV0dXJuIFtdXG4gIH1cbiAgXG4gIHZhciBkID0gcG9pbnRzWzBdLmxlbmd0aFxuICBpZihkIDwgMSkge1xuICAgIHJldHVybiBbXVxuICB9XG5cbiAgLy9TcGVjaWFsIGNhc2U6ICBGb3IgMUQgd2UgY2FuIGp1c3Qgc29ydCB0aGUgcG9pbnRzXG4gIGlmKGQgPT09IDEpIHtcbiAgICByZXR1cm4gdHJpYW5ndWxhdGUxRChuLCBwb2ludHMsIGluY2x1ZGVQb2ludEF0SW5maW5pdHkpXG4gIH1cbiAgXG4gIC8vTGlmdCBwb2ludHMsIHNvcnRcbiAgdmFyIGxpZnRlZCA9IG5ldyBBcnJheShuKVxuICB2YXIgdXBwZXIgPSAxLjBcbiAgZm9yKHZhciBpPTA7IGk8bjsgKytpKSB7XG4gICAgdmFyIHAgPSBwb2ludHNbaV1cbiAgICB2YXIgeCA9IG5ldyBBcnJheShkKzEpXG4gICAgdmFyIGwgPSAwLjBcbiAgICBmb3IodmFyIGo9MDsgajxkOyArK2opIHtcbiAgICAgIHZhciB2ID0gcFtqXVxuICAgICAgeFtqXSA9IHZcbiAgICAgIGwgKz0gdiAqIHZcbiAgICB9XG4gICAgeFtkXSA9IGxcbiAgICBsaWZ0ZWRbaV0gPSBuZXcgTGlmdGVkUG9pbnQoeCwgaSlcbiAgICB1cHBlciA9IE1hdGgubWF4KGwsIHVwcGVyKVxuICB9XG4gIHVuaXEobGlmdGVkLCBjb21wYXJlTGlmdGVkKVxuICBcbiAgLy9Eb3VibGUgcG9pbnRzXG4gIG4gPSBsaWZ0ZWQubGVuZ3RoXG5cbiAgLy9DcmVhdGUgbmV3IGxpc3Qgb2YgcG9pbnRzXG4gIHZhciBkcG9pbnRzID0gbmV3IEFycmF5KG4gKyBkICsgMSlcbiAgdmFyIGRpbmRleCA9IG5ldyBBcnJheShuICsgZCArIDEpXG5cbiAgLy9BZGQgc3RlaW5lciBwb2ludHMgYXQgdG9wXG4gIHZhciB1ID0gKGQrMSkgKiAoZCsxKSAqIHVwcGVyXG4gIHZhciB5ID0gbmV3IEFycmF5KGQrMSlcbiAgZm9yKHZhciBpPTA7IGk8PWQ7ICsraSkge1xuICAgIHlbaV0gPSAwLjBcbiAgfVxuICB5W2RdID0gdVxuXG4gIGRwb2ludHNbMF0gPSB5LnNsaWNlKClcbiAgZGluZGV4WzBdID0gLTFcblxuICBmb3IodmFyIGk9MDsgaTw9ZDsgKytpKSB7XG4gICAgdmFyIHggPSB5LnNsaWNlKClcbiAgICB4W2ldID0gMVxuICAgIGRwb2ludHNbaSsxXSA9IHhcbiAgICBkaW5kZXhbaSsxXSA9IC0xXG4gIH1cblxuICAvL0NvcHkgcmVzdCBvZiB0aGUgcG9pbnRzIG92ZXJcbiAgZm9yKHZhciBpPTA7IGk8bjsgKytpKSB7XG4gICAgdmFyIGggPSBsaWZ0ZWRbaV1cbiAgICBkcG9pbnRzW2kgKyBkICsgMV0gPSBoLnBvaW50XG4gICAgZGluZGV4W2kgKyBkICsgMV0gPSAgaC5pbmRleFxuICB9XG5cbiAgLy9Db25zdHJ1Y3QgY29udmV4IGh1bGxcbiAgdmFyIGh1bGwgPSBjaChkcG9pbnRzLCBmYWxzZSlcbiAgaWYoaW5jbHVkZVBvaW50QXRJbmZpbml0eSkge1xuICAgIGh1bGwgPSBodWxsLmZpbHRlcihmdW5jdGlvbihjZWxsKSB7XG4gICAgICB2YXIgY291bnQgPSAwXG4gICAgICBmb3IodmFyIGo9MDsgajw9ZDsgKytqKSB7XG4gICAgICAgIHZhciB2ID0gZGluZGV4W2NlbGxbal1dXG4gICAgICAgIGlmKHYgPCAwKSB7XG4gICAgICAgICAgaWYoKytjb3VudCA+PSAyKSB7XG4gICAgICAgICAgICByZXR1cm4gZmFsc2VcbiAgICAgICAgICB9XG4gICAgICAgIH1cbiAgICAgICAgY2VsbFtqXSA9IHZcbiAgICAgIH1cbiAgICAgIHJldHVybiB0cnVlXG4gICAgfSlcbiAgfSBlbHNlIHtcbiAgICBodWxsID0gaHVsbC5maWx0ZXIoZnVuY3Rpb24oY2VsbCkge1xuICAgICAgZm9yKHZhciBpPTA7IGk8PWQ7ICsraSkge1xuICAgICAgICB2YXIgdiA9IGRpbmRleFtjZWxsW2ldXVxuICAgICAgICBpZih2IDwgMCkge1xuICAgICAgICAgIHJldHVybiBmYWxzZVxuICAgICAgICB9XG4gICAgICAgIGNlbGxbaV0gPSB2XG4gICAgICB9XG4gICAgICByZXR1cm4gdHJ1ZVxuICAgIH0pXG4gIH1cblxuICBpZihkICYgMSkge1xuICAgIGZvcih2YXIgaT0wOyBpPGh1bGwubGVuZ3RoOyArK2kpIHtcbiAgICAgIHZhciBoID0gaHVsbFtpXVxuICAgICAgdmFyIHggPSBoWzBdXG4gICAgICBoWzBdID0gaFsxXVxuICAgICAgaFsxXSA9IHhcbiAgICB9XG4gIH1cblxuICByZXR1cm4gaHVsbFxufSIsImlmKHR5cGVvZiB3aW5kb3cucGVyZm9ybWFuY2UgPT09IFwib2JqZWN0XCIpIHtcbiAgaWYod2luZG93LnBlcmZvcm1hbmNlLm5vdykge1xuICAgIG1vZHVsZS5leHBvcnRzID0gZnVuY3Rpb24oKSB7IHJldHVybiB3aW5kb3cucGVyZm9ybWFuY2Uubm93KCkgfVxuICB9IGVsc2UgaWYod2luZG93LnBlcmZvcm1hbmNlLndlYmtpdE5vdykge1xuICAgIG1vZHVsZS5leHBvcnRzID0gZnVuY3Rpb24oKSB7IHJldHVybiB3aW5kb3cucGVyZm9ybWFuY2Uud2Via2l0Tm93KCkgfVxuICB9XG59IGVsc2UgaWYoRGF0ZS5ub3cpIHtcbiAgbW9kdWxlLmV4cG9ydHMgPSBEYXRlLm5vd1xufSBlbHNlIHtcbiAgbW9kdWxlLmV4cG9ydHMgPSBmdW5jdGlvbigpIHsgcmV0dXJuIChuZXcgRGF0ZSgpKS5nZXRUaW1lKCkgfVxufVxuIiwiLy9BZGFwdGVkIGZyb20gaGVyZTogaHR0cHM6Ly9kZXZlbG9wZXIubW96aWxsYS5vcmcvZW4tVVMvZG9jcy9XZWIvUmVmZXJlbmNlL0V2ZW50cy93aGVlbD9yZWRpcmVjdGxvY2FsZT1lbi1VUyZyZWRpcmVjdHNsdWc9RE9NJTJGTW96aWxsYV9ldmVudF9yZWZlcmVuY2UlMkZ3aGVlbFxuXG52YXIgcHJlZml4ID0gXCJcIiwgX2FkZEV2ZW50TGlzdGVuZXIsIG9ud2hlZWwsIHN1cHBvcnQ7XG5cbi8vIGRldGVjdCBldmVudCBtb2RlbFxuaWYgKCB3aW5kb3cuYWRkRXZlbnRMaXN0ZW5lciApIHtcbiAgX2FkZEV2ZW50TGlzdGVuZXIgPSBcImFkZEV2ZW50TGlzdGVuZXJcIjtcbn0gZWxzZSB7XG4gIF9hZGRFdmVudExpc3RlbmVyID0gXCJhdHRhY2hFdmVudFwiO1xuICBwcmVmaXggPSBcIm9uXCI7XG59XG5cbi8vIGRldGVjdCBhdmFpbGFibGUgd2hlZWwgZXZlbnRcbnN1cHBvcnQgPSBcIm9ud2hlZWxcIiBpbiBkb2N1bWVudC5jcmVhdGVFbGVtZW50KFwiZGl2XCIpID8gXCJ3aGVlbFwiIDogLy8gTW9kZXJuIGJyb3dzZXJzIHN1cHBvcnQgXCJ3aGVlbFwiXG4gICAgICAgICAgZG9jdW1lbnQub25tb3VzZXdoZWVsICE9PSB1bmRlZmluZWQgPyBcIm1vdXNld2hlZWxcIiA6IC8vIFdlYmtpdCBhbmQgSUUgc3VwcG9ydCBhdCBsZWFzdCBcIm1vdXNld2hlZWxcIlxuICAgICAgICAgIFwiRE9NTW91c2VTY3JvbGxcIjsgLy8gbGV0J3MgYXNzdW1lIHRoYXQgcmVtYWluaW5nIGJyb3dzZXJzIGFyZSBvbGRlciBGaXJlZm94XG5cbmZ1bmN0aW9uIF9hZGRXaGVlbExpc3RlbmVyKCBlbGVtLCBldmVudE5hbWUsIGNhbGxiYWNrLCB1c2VDYXB0dXJlICkge1xuICBlbGVtWyBfYWRkRXZlbnRMaXN0ZW5lciBdKCBwcmVmaXggKyBldmVudE5hbWUsIHN1cHBvcnQgPT0gXCJ3aGVlbFwiID8gY2FsbGJhY2sgOiBmdW5jdGlvbiggb3JpZ2luYWxFdmVudCApIHtcbiAgICAhb3JpZ2luYWxFdmVudCAmJiAoIG9yaWdpbmFsRXZlbnQgPSB3aW5kb3cuZXZlbnQgKTtcblxuICAgIC8vIGNyZWF0ZSBhIG5vcm1hbGl6ZWQgZXZlbnQgb2JqZWN0XG4gICAgdmFyIGV2ZW50ID0ge1xuICAgICAgLy8ga2VlcCBhIHJlZiB0byB0aGUgb3JpZ2luYWwgZXZlbnQgb2JqZWN0XG4gICAgICBvcmlnaW5hbEV2ZW50OiBvcmlnaW5hbEV2ZW50LFxuICAgICAgdGFyZ2V0OiBvcmlnaW5hbEV2ZW50LnRhcmdldCB8fCBvcmlnaW5hbEV2ZW50LnNyY0VsZW1lbnQsXG4gICAgICB0eXBlOiBcIndoZWVsXCIsXG4gICAgICBkZWx0YU1vZGU6IG9yaWdpbmFsRXZlbnQudHlwZSA9PSBcIk1vek1vdXNlUGl4ZWxTY3JvbGxcIiA/IDAgOiAxLFxuICAgICAgZGVsdGFYOiAwLFxuICAgICAgZGVsYXRaOiAwLFxuICAgICAgcHJldmVudERlZmF1bHQ6IGZ1bmN0aW9uKCkge1xuICAgICAgICBvcmlnaW5hbEV2ZW50LnByZXZlbnREZWZhdWx0ID9cbiAgICAgICAgICBvcmlnaW5hbEV2ZW50LnByZXZlbnREZWZhdWx0KCkgOlxuICAgICAgICAgIG9yaWdpbmFsRXZlbnQucmV0dXJuVmFsdWUgPSBmYWxzZTtcbiAgICAgIH1cbiAgICB9O1xuICAgIFxuICAgIC8vIGNhbGN1bGF0ZSBkZWx0YVkgKGFuZCBkZWx0YVgpIGFjY29yZGluZyB0byB0aGUgZXZlbnRcbiAgICBpZiAoIHN1cHBvcnQgPT0gXCJtb3VzZXdoZWVsXCIgKSB7XG4gICAgICBldmVudC5kZWx0YVkgPSAtIDEvNDAgKiBvcmlnaW5hbEV2ZW50LndoZWVsRGVsdGE7XG4gICAgICAvLyBXZWJraXQgYWxzbyBzdXBwb3J0IHdoZWVsRGVsdGFYXG4gICAgICBvcmlnaW5hbEV2ZW50LndoZWVsRGVsdGFYICYmICggZXZlbnQuZGVsdGFYID0gLSAxLzQwICogb3JpZ2luYWxFdmVudC53aGVlbERlbHRhWCApO1xuICAgIH0gZWxzZSB7XG4gICAgICBldmVudC5kZWx0YVkgPSBvcmlnaW5hbEV2ZW50LmRldGFpbDtcbiAgICB9XG5cbiAgICAvLyBpdCdzIHRpbWUgdG8gZmlyZSB0aGUgY2FsbGJhY2tcbiAgICByZXR1cm4gY2FsbGJhY2soIGV2ZW50ICk7XG4gIH0sIHVzZUNhcHR1cmUgfHwgZmFsc2UgKTtcbn1cblxubW9kdWxlLmV4cG9ydHMgPSBmdW5jdGlvbiggZWxlbSwgY2FsbGJhY2ssIHVzZUNhcHR1cmUgKSB7XG4gIF9hZGRXaGVlbExpc3RlbmVyKCBlbGVtLCBzdXBwb3J0LCBjYWxsYmFjaywgdXNlQ2FwdHVyZSApO1xuXG4gIC8vIGhhbmRsZSBNb3pNb3VzZVBpeGVsU2Nyb2xsIGluIG9sZGVyIEZpcmVmb3hcbiAgaWYoIHN1cHBvcnQgPT0gXCJET01Nb3VzZVNjcm9sbFwiICkge1xuICAgIF9hZGRXaGVlbExpc3RlbmVyKCBlbGVtLCBcIk1vek1vdXNlUGl4ZWxTY3JvbGxcIiwgY2FsbGJhY2ssIHVzZUNhcHR1cmUgKTtcbiAgfVxufTsiLCIvLyBodHRwOi8vcGF1bGlyaXNoLmNvbS8yMDExL3JlcXVlc3RhbmltYXRpb25mcmFtZS1mb3Itc21hcnQtYW5pbWF0aW5nL1xuLy8gaHR0cDovL215Lm9wZXJhLmNvbS9lbW9sbGVyL2Jsb2cvMjAxMS8xMi8yMC9yZXF1ZXN0YW5pbWF0aW9uZnJhbWUtZm9yLXNtYXJ0LWVyLWFuaW1hdGluZ1xuIFxuLy8gcmVxdWVzdEFuaW1hdGlvbkZyYW1lIHBvbHlmaWxsIGJ5IEVyaWsgTcO2bGxlci4gZml4ZXMgZnJvbSBQYXVsIElyaXNoIGFuZCBUaW5vIFppamRlbFxuIFxuLy8gTUlUIGxpY2Vuc2VcbnZhciBsYXN0VGltZSA9IDA7XG52YXIgdmVuZG9ycyA9IFsnbXMnLCAnbW96JywgJ3dlYmtpdCcsICdvJ107XG5mb3IodmFyIHggPSAwOyB4IDwgdmVuZG9ycy5sZW5ndGggJiYgIXdpbmRvdy5yZXF1ZXN0QW5pbWF0aW9uRnJhbWU7ICsreCkge1xuICAgIHdpbmRvdy5yZXF1ZXN0QW5pbWF0aW9uRnJhbWUgPSB3aW5kb3dbdmVuZG9yc1t4XSsnUmVxdWVzdEFuaW1hdGlvbkZyYW1lJ107XG4gICAgd2luZG93LmNhbmNlbEFuaW1hdGlvbkZyYW1lID0gd2luZG93W3ZlbmRvcnNbeF0rJ0NhbmNlbEFuaW1hdGlvbkZyYW1lJ10gXG4gICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgfHwgd2luZG93W3ZlbmRvcnNbeF0rJ0NhbmNlbFJlcXVlc3RBbmltYXRpb25GcmFtZSddO1xufVxuXG5pZiAoIXdpbmRvdy5yZXF1ZXN0QW5pbWF0aW9uRnJhbWUpXG4gICAgd2luZG93LnJlcXVlc3RBbmltYXRpb25GcmFtZSA9IGZ1bmN0aW9uKGNhbGxiYWNrLCBlbGVtZW50KSB7XG4gICAgICAgIHZhciBjdXJyVGltZSA9IG5ldyBEYXRlKCkuZ2V0VGltZSgpO1xuICAgICAgICB2YXIgdGltZVRvQ2FsbCA9IE1hdGgubWF4KDAsIDE2IC0gKGN1cnJUaW1lIC0gbGFzdFRpbWUpKTtcbiAgICAgICAgdmFyIGlkID0gd2luZG93LnNldFRpbWVvdXQoZnVuY3Rpb24oKSB7IGNhbGxiYWNrKGN1cnJUaW1lICsgdGltZVRvQ2FsbCk7IH0sIFxuICAgICAgICAgIHRpbWVUb0NhbGwpO1xuICAgICAgICBsYXN0VGltZSA9IGN1cnJUaW1lICsgdGltZVRvQ2FsbDtcbiAgICAgICAgcmV0dXJuIGlkO1xuICAgIH07XG5cbmlmICghd2luZG93LmNhbmNlbEFuaW1hdGlvbkZyYW1lKVxuICAgIHdpbmRvdy5jYW5jZWxBbmltYXRpb25GcmFtZSA9IGZ1bmN0aW9uKGlkKSB7XG4gICAgICAgIGNsZWFyVGltZW91dChpZCk7XG4gICAgfTtcbiIsIlwidXNlIHN0cmljdFwiXG5cbmZ1bmN0aW9uIGNvbXBpbGVTZWFyY2goZnVuY05hbWUsIHByZWRpY2F0ZSwgcmV2ZXJzZWQsIGV4dHJhQXJncywgdXNlTmRhcnJheSwgZWFybHlPdXQpIHtcbiAgdmFyIGNvZGUgPSBbXG4gICAgXCJmdW5jdGlvbiBcIiwgZnVuY05hbWUsIFwiKGEsbCxoLFwiLCBleHRyYUFyZ3Muam9pbihcIixcIiksICBcIil7XCIsXG5lYXJseU91dCA/IFwiXCIgOiBcInZhciBpPVwiLCAocmV2ZXJzZWQgPyBcImwtMVwiIDogXCJoKzFcIiksXG5cIjt3aGlsZShsPD1oKXtcXFxudmFyIG09KGwraCk+Pj4xLHg9YVwiLCB1c2VOZGFycmF5ID8gXCIuZ2V0KG0pXCIgOiBcIlttXVwiXVxuICBpZihlYXJseU91dCkge1xuICAgIGlmKHByZWRpY2F0ZS5pbmRleE9mKFwiY1wiKSA8IDApIHtcbiAgICAgIGNvZGUucHVzaChcIjtpZih4PT09eSl7cmV0dXJuIG19ZWxzZSBpZih4PD15KXtcIilcbiAgICB9IGVsc2Uge1xuICAgICAgY29kZS5wdXNoKFwiO3ZhciBwPWMoeCx5KTtpZihwPT09MCl7cmV0dXJuIG19ZWxzZSBpZihwPD0wKXtcIilcbiAgICB9XG4gIH0gZWxzZSB7XG4gICAgY29kZS5wdXNoKFwiO2lmKFwiLCBwcmVkaWNhdGUsIFwiKXtpPW07XCIpXG4gIH1cbiAgaWYocmV2ZXJzZWQpIHtcbiAgICBjb2RlLnB1c2goXCJsPW0rMX1lbHNle2g9bS0xfVwiKVxuICB9IGVsc2Uge1xuICAgIGNvZGUucHVzaChcImg9bS0xfWVsc2V7bD1tKzF9XCIpXG4gIH1cbiAgY29kZS5wdXNoKFwifVwiKVxuICBpZihlYXJseU91dCkge1xuICAgIGNvZGUucHVzaChcInJldHVybiAtMX07XCIpXG4gIH0gZWxzZSB7XG4gICAgY29kZS5wdXNoKFwicmV0dXJuIGl9O1wiKVxuICB9XG4gIHJldHVybiBjb2RlLmpvaW4oXCJcIilcbn1cblxuZnVuY3Rpb24gY29tcGlsZUJvdW5kc1NlYXJjaChwcmVkaWNhdGUsIHJldmVyc2VkLCBzdWZmaXgsIGVhcmx5T3V0KSB7XG4gIHZhciByZXN1bHQgPSBuZXcgRnVuY3Rpb24oW1xuICBjb21waWxlU2VhcmNoKFwiQVwiLCBcInhcIiArIHByZWRpY2F0ZSArIFwieVwiLCByZXZlcnNlZCwgW1wieVwiXSwgZmFsc2UsIGVhcmx5T3V0KSxcbiAgY29tcGlsZVNlYXJjaChcIkJcIiwgXCJ4XCIgKyBwcmVkaWNhdGUgKyBcInlcIiwgcmV2ZXJzZWQsIFtcInlcIl0sIHRydWUsIGVhcmx5T3V0KSxcbiAgY29tcGlsZVNlYXJjaChcIlBcIiwgXCJjKHgseSlcIiArIHByZWRpY2F0ZSArIFwiMFwiLCByZXZlcnNlZCwgW1wieVwiLCBcImNcIl0sIGZhbHNlLCBlYXJseU91dCksXG4gIGNvbXBpbGVTZWFyY2goXCJRXCIsIFwiYyh4LHkpXCIgKyBwcmVkaWNhdGUgKyBcIjBcIiwgcmV2ZXJzZWQsIFtcInlcIiwgXCJjXCJdLCB0cnVlLCBlYXJseU91dCksXG5cImZ1bmN0aW9uIGRpc3BhdGNoQnNlYXJjaFwiLCBzdWZmaXgsIFwiKGEseSxjLGwsaCl7XFxcbmlmKGEuc2hhcGUpe1xcXG5pZih0eXBlb2YoYyk9PT0nZnVuY3Rpb24nKXtcXFxucmV0dXJuIFEoYSwobD09PXVuZGVmaW5lZCk/MDpsfDAsKGg9PT11bmRlZmluZWQpP2Euc2hhcGVbMF0tMTpofDAseSxjKVxcXG59ZWxzZXtcXFxucmV0dXJuIEIoYSwoYz09PXVuZGVmaW5lZCk/MDpjfDAsKGw9PT11bmRlZmluZWQpP2Euc2hhcGVbMF0tMTpsfDAseSlcXFxufX1lbHNle1xcXG5pZih0eXBlb2YoYyk9PT0nZnVuY3Rpb24nKXtcXFxucmV0dXJuIFAoYSwobD09PXVuZGVmaW5lZCk/MDpsfDAsKGg9PT11bmRlZmluZWQpP2EubGVuZ3RoLTE6aHwwLHksYylcXFxufWVsc2V7XFxcbnJldHVybiBBKGEsKGM9PT11bmRlZmluZWQpPzA6Y3wwLChsPT09dW5kZWZpbmVkKT9hLmxlbmd0aC0xOmx8MCx5KVxcXG59fX1cXFxucmV0dXJuIGRpc3BhdGNoQnNlYXJjaFwiLCBzdWZmaXhdLmpvaW4oXCJcIikpXG4gIHJldHVybiByZXN1bHQoKVxufVxuXG5tb2R1bGUuZXhwb3J0cyA9IHtcbiAgZ2U6IGNvbXBpbGVCb3VuZHNTZWFyY2goXCI+PVwiLCBmYWxzZSwgXCJHRVwiKSxcbiAgZ3Q6IGNvbXBpbGVCb3VuZHNTZWFyY2goXCI+XCIsIGZhbHNlLCBcIkdUXCIpLFxuICBsdDogY29tcGlsZUJvdW5kc1NlYXJjaChcIjxcIiwgdHJ1ZSwgXCJMVFwiKSxcbiAgbGU6IGNvbXBpbGVCb3VuZHNTZWFyY2goXCI8PVwiLCB0cnVlLCBcIkxFXCIpLFxuICBlcTogY29tcGlsZUJvdW5kc1NlYXJjaChcIi1cIiwgdHJ1ZSwgXCJFUVwiLCB0cnVlKVxufVxuIiwiLyohXG4gICogZG9tcmVhZHkgKGMpIER1c3RpbiBEaWF6IDIwMTQgLSBMaWNlbnNlIE1JVFxuICAqL1xuIWZ1bmN0aW9uIChuYW1lLCBkZWZpbml0aW9uKSB7XG5cbiAgaWYgKHR5cGVvZiBtb2R1bGUgIT0gJ3VuZGVmaW5lZCcpIG1vZHVsZS5leHBvcnRzID0gZGVmaW5pdGlvbigpXG4gIGVsc2UgaWYgKHR5cGVvZiBkZWZpbmUgPT0gJ2Z1bmN0aW9uJyAmJiB0eXBlb2YgZGVmaW5lLmFtZCA9PSAnb2JqZWN0JykgZGVmaW5lKGRlZmluaXRpb24pXG4gIGVsc2UgdGhpc1tuYW1lXSA9IGRlZmluaXRpb24oKVxuXG59KCdkb21yZWFkeScsIGZ1bmN0aW9uICgpIHtcblxuICB2YXIgZm5zID0gW10sIGxpc3RlbmVyXG4gICAgLCBkb2MgPSBkb2N1bWVudFxuICAgICwgZG9tQ29udGVudExvYWRlZCA9ICdET01Db250ZW50TG9hZGVkJ1xuICAgICwgbG9hZGVkID0gL15sb2FkZWR8XmMvLnRlc3QoZG9jLnJlYWR5U3RhdGUpXG5cbiAgaWYgKCFsb2FkZWQpXG4gIGRvYy5hZGRFdmVudExpc3RlbmVyKGRvbUNvbnRlbnRMb2FkZWQsIGxpc3RlbmVyID0gZnVuY3Rpb24gKCkge1xuICAgIGRvYy5yZW1vdmVFdmVudExpc3RlbmVyKGRvbUNvbnRlbnRMb2FkZWQsIGxpc3RlbmVyKVxuICAgIGxvYWRlZCA9IDFcbiAgICB3aGlsZSAobGlzdGVuZXIgPSBmbnMuc2hpZnQoKSkgbGlzdGVuZXIoKVxuICB9KVxuXG4gIHJldHVybiBmdW5jdGlvbiAoZm4pIHtcbiAgICBsb2FkZWQgPyBmbigpIDogZm5zLnB1c2goZm4pXG4gIH1cblxufSk7XG4iLCJcInVzZSBzdHJpY3RcIlxuXG5mdW5jdGlvbiBpbnZlcnQoaGFzaCkge1xuICB2YXIgcmVzdWx0ID0ge31cbiAgZm9yKHZhciBpIGluIGhhc2gpIHtcbiAgICBpZihoYXNoLmhhc093blByb3BlcnR5KGkpKSB7XG4gICAgICByZXN1bHRbaGFzaFtpXV0gPSBpXG4gICAgfVxuICB9XG4gIHJldHVybiByZXN1bHRcbn1cblxubW9kdWxlLmV4cG9ydHMgPSBpbnZlcnQiLCJcInVzZSBzdHJpY3RcIlxuXG5mdW5jdGlvbiBpb3RhKG4pIHtcbiAgdmFyIHJlc3VsdCA9IG5ldyBBcnJheShuKVxuICBmb3IodmFyIGk9MDsgaTxuOyArK2kpIHtcbiAgICByZXN1bHRbaV0gPSBpXG4gIH1cbiAgcmV0dXJuIHJlc3VsdFxufVxuXG5tb2R1bGUuZXhwb3J0cyA9IGlvdGEiLCJcInVzZSBzdHJpY3RcIlxuXG5mdW5jdGlvbiB1bmlxdWVfcHJlZChsaXN0LCBjb21wYXJlKSB7XG4gIHZhciBwdHIgPSAxXG4gICAgLCBsZW4gPSBsaXN0Lmxlbmd0aFxuICAgICwgYT1saXN0WzBdLCBiPWxpc3RbMF1cbiAgZm9yKHZhciBpPTE7IGk8bGVuOyArK2kpIHtcbiAgICBiID0gYVxuICAgIGEgPSBsaXN0W2ldXG4gICAgaWYoY29tcGFyZShhLCBiKSkge1xuICAgICAgaWYoaSA9PT0gcHRyKSB7XG4gICAgICAgIHB0cisrXG4gICAgICAgIGNvbnRpbnVlXG4gICAgICB9XG4gICAgICBsaXN0W3B0cisrXSA9IGFcbiAgICB9XG4gIH1cbiAgbGlzdC5sZW5ndGggPSBwdHJcbiAgcmV0dXJuIGxpc3Rcbn1cblxuZnVuY3Rpb24gdW5pcXVlX2VxKGxpc3QpIHtcbiAgdmFyIHB0ciA9IDFcbiAgICAsIGxlbiA9IGxpc3QubGVuZ3RoXG4gICAgLCBhPWxpc3RbMF0sIGIgPSBsaXN0WzBdXG4gIGZvcih2YXIgaT0xOyBpPGxlbjsgKytpLCBiPWEpIHtcbiAgICBiID0gYVxuICAgIGEgPSBsaXN0W2ldXG4gICAgaWYoYSAhPT0gYikge1xuICAgICAgaWYoaSA9PT0gcHRyKSB7XG4gICAgICAgIHB0cisrXG4gICAgICAgIGNvbnRpbnVlXG4gICAgICB9XG4gICAgICBsaXN0W3B0cisrXSA9IGFcbiAgICB9XG4gIH1cbiAgbGlzdC5sZW5ndGggPSBwdHJcbiAgcmV0dXJuIGxpc3Rcbn1cblxuZnVuY3Rpb24gdW5pcXVlKGxpc3QsIGNvbXBhcmUsIHNvcnRlZCkge1xuICBpZihsaXN0Lmxlbmd0aCA9PT0gMCkge1xuICAgIHJldHVybiBsaXN0XG4gIH1cbiAgaWYoY29tcGFyZSkge1xuICAgIGlmKCFzb3J0ZWQpIHtcbiAgICAgIGxpc3Quc29ydChjb21wYXJlKVxuICAgIH1cbiAgICByZXR1cm4gdW5pcXVlX3ByZWQobGlzdCwgY29tcGFyZSlcbiAgfVxuICBpZighc29ydGVkKSB7XG4gICAgbGlzdC5zb3J0KClcbiAgfVxuICByZXR1cm4gdW5pcXVlX2VxKGxpc3QpXG59XG5cbm1vZHVsZS5leHBvcnRzID0gdW5pcXVlXG4iLCJ2YXIgdWEgPSB0eXBlb2Ygd2luZG93ICE9PSAndW5kZWZpbmVkJyA/IHdpbmRvdy5uYXZpZ2F0b3IudXNlckFnZW50IDogJydcbiAgLCBpc09TWCA9IC9PUyBYLy50ZXN0KHVhKVxuICAsIGlzT3BlcmEgPSAvT3BlcmEvLnRlc3QodWEpXG4gICwgbWF5YmVGaXJlZm94ID0gIS9saWtlIEdlY2tvLy50ZXN0KHVhKSAmJiAhaXNPcGVyYVxuXG52YXIgaSwgb3V0cHV0ID0gbW9kdWxlLmV4cG9ydHMgPSB7XG4gIDA6ICBpc09TWCA/ICc8bWVudT4nIDogJzxVTks+J1xuLCAxOiAgJzxtb3VzZSAxPidcbiwgMjogICc8bW91c2UgMj4nXG4sIDM6ICAnPGJyZWFrPidcbiwgNDogICc8bW91c2UgMz4nXG4sIDU6ICAnPG1vdXNlIDQ+J1xuLCA2OiAgJzxtb3VzZSA1PidcbiwgODogICc8YmFja3NwYWNlPidcbiwgOTogICc8dGFiPidcbiwgMTI6ICc8Y2xlYXI+J1xuLCAxMzogJzxlbnRlcj4nXG4sIDE2OiAnPHNoaWZ0PidcbiwgMTc6ICc8Y29udHJvbD4nXG4sIDE4OiAnPGFsdD4nXG4sIDE5OiAnPHBhdXNlPidcbiwgMjA6ICc8Y2Fwcy1sb2NrPidcbiwgMjE6ICc8aW1lLWhhbmd1bD4nXG4sIDIzOiAnPGltZS1qdW5qYT4nXG4sIDI0OiAnPGltZS1maW5hbD4nXG4sIDI1OiAnPGltZS1rYW5qaT4nXG4sIDI3OiAnPGVzY2FwZT4nXG4sIDI4OiAnPGltZS1jb252ZXJ0PidcbiwgMjk6ICc8aW1lLW5vbmNvbnZlcnQ+J1xuLCAzMDogJzxpbWUtYWNjZXB0PidcbiwgMzE6ICc8aW1lLW1vZGUtY2hhbmdlPidcbiwgMjc6ICc8ZXNjYXBlPidcbiwgMzI6ICc8c3BhY2U+J1xuLCAzMzogJzxwYWdlLXVwPidcbiwgMzQ6ICc8cGFnZS1kb3duPidcbiwgMzU6ICc8ZW5kPidcbiwgMzY6ICc8aG9tZT4nXG4sIDM3OiAnPGxlZnQ+J1xuLCAzODogJzx1cD4nXG4sIDM5OiAnPHJpZ2h0PidcbiwgNDA6ICc8ZG93bj4nXG4sIDQxOiAnPHNlbGVjdD4nXG4sIDQyOiAnPHByaW50PidcbiwgNDM6ICc8ZXhlY3V0ZT4nXG4sIDQ0OiAnPHNuYXBzaG90PidcbiwgNDU6ICc8aW5zZXJ0PidcbiwgNDY6ICc8ZGVsZXRlPidcbiwgNDc6ICc8aGVscD4nXG4sIDkxOiAnPG1ldGE+JyAgLy8gbWV0YS1sZWZ0IC0tIG5vIG9uZSBoYW5kbGVzIGxlZnQgYW5kIHJpZ2h0IHByb3Blcmx5LCBzbyB3ZSBjb2VyY2UgaW50byBvbmUuXG4sIDkyOiAnPG1ldGE+JyAgLy8gbWV0YS1yaWdodFxuLCA5MzogaXNPU1ggPyAnPG1ldGE+JyA6ICc8bWVudT4nICAgICAgLy8gY2hyb21lLG9wZXJhLHNhZmFyaSBhbGwgcmVwb3J0IHRoaXMgZm9yIG1ldGEtcmlnaHQgKG9zeCBtYnApLlxuLCA5NTogJzxzbGVlcD4nXG4sIDEwNjogJzxudW0tKj4nXG4sIDEwNzogJzxudW0tKz4nXG4sIDEwODogJzxudW0tZW50ZXI+J1xuLCAxMDk6ICc8bnVtLS0+J1xuLCAxMTA6ICc8bnVtLS4+J1xuLCAxMTE6ICc8bnVtLS8+J1xuLCAxNDQ6ICc8bnVtLWxvY2s+J1xuLCAxNDU6ICc8c2Nyb2xsLWxvY2s+J1xuLCAxNjA6ICc8c2hpZnQtbGVmdD4nXG4sIDE2MTogJzxzaGlmdC1yaWdodD4nXG4sIDE2MjogJzxjb250cm9sLWxlZnQ+J1xuLCAxNjM6ICc8Y29udHJvbC1yaWdodD4nXG4sIDE2NDogJzxhbHQtbGVmdD4nXG4sIDE2NTogJzxhbHQtcmlnaHQ+J1xuLCAxNjY6ICc8YnJvd3Nlci1iYWNrPidcbiwgMTY3OiAnPGJyb3dzZXItZm9yd2FyZD4nXG4sIDE2ODogJzxicm93c2VyLXJlZnJlc2g+J1xuLCAxNjk6ICc8YnJvd3Nlci1zdG9wPidcbiwgMTcwOiAnPGJyb3dzZXItc2VhcmNoPidcbiwgMTcxOiAnPGJyb3dzZXItZmF2b3JpdGVzPidcbiwgMTcyOiAnPGJyb3dzZXItaG9tZT4nXG5cbiAgLy8gZmYvb3N4IHJlcG9ydHMgJzx2b2x1bWUtbXV0ZT4nIGZvciAnLSdcbiwgMTczOiBpc09TWCAmJiBtYXliZUZpcmVmb3ggPyAnLScgOiAnPHZvbHVtZS1tdXRlPidcbiwgMTc0OiAnPHZvbHVtZS1kb3duPidcbiwgMTc1OiAnPHZvbHVtZS11cD4nXG4sIDE3NjogJzxuZXh0LXRyYWNrPidcbiwgMTc3OiAnPHByZXYtdHJhY2s+J1xuLCAxNzg6ICc8c3RvcD4nXG4sIDE3OTogJzxwbGF5LXBhdXNlPidcbiwgMTgwOiAnPGxhdW5jaC1tYWlsPidcbiwgMTgxOiAnPGxhdW5jaC1tZWRpYS1zZWxlY3Q+J1xuLCAxODI6ICc8bGF1bmNoLWFwcCAxPidcbiwgMTgzOiAnPGxhdW5jaC1hcHAgMj4nXG4sIDE4NjogJzsnXG4sIDE4NzogJz0nXG4sIDE4ODogJywnXG4sIDE4OTogJy0nXG4sIDE5MDogJy4nXG4sIDE5MTogJy8nXG4sIDE5MjogJ2AnXG4sIDIxOTogJ1snXG4sIDIyMDogJ1xcXFwnXG4sIDIyMTogJ10nXG4sIDIyMjogXCInXCJcbiwgMjIzOiAnPG1ldGE+J1xuLCAyMjQ6ICc8bWV0YT4nICAgICAgIC8vIGZpcmVmb3ggcmVwb3J0cyBtZXRhIGhlcmUuXG4sIDIyNjogJzxhbHQtZ3I+J1xuLCAyMjk6ICc8aW1lLXByb2Nlc3M+J1xuLCAyMzE6IGlzT3BlcmEgPyAnYCcgOiAnPHVuaWNvZGU+J1xuLCAyNDY6ICc8YXR0ZW50aW9uPidcbiwgMjQ3OiAnPGNyc2VsPidcbiwgMjQ4OiAnPGV4c2VsPidcbiwgMjQ5OiAnPGVyYXNlLWVvZj4nXG4sIDI1MDogJzxwbGF5PidcbiwgMjUxOiAnPHpvb20+J1xuLCAyNTI6ICc8bm8tbmFtZT4nXG4sIDI1MzogJzxwYS0xPidcbiwgMjU0OiAnPGNsZWFyPidcbn1cblxuZm9yKGkgPSA1ODsgaSA8IDY1OyArK2kpIHtcbiAgb3V0cHV0W2ldID0gU3RyaW5nLmZyb21DaGFyQ29kZShpKVxufVxuXG4vLyAwLTlcbmZvcihpID0gNDg7IGkgPCA1ODsgKytpKSB7XG4gIG91dHB1dFtpXSA9IChpIC0gNDgpKycnXG59XG5cbi8vIEEtWlxuZm9yKGkgPSA2NTsgaSA8IDkxOyArK2kpIHtcbiAgb3V0cHV0W2ldID0gU3RyaW5nLmZyb21DaGFyQ29kZShpKVxufVxuXG4vLyBudW0wLTlcbmZvcihpID0gOTY7IGkgPCAxMDY7ICsraSkge1xuICBvdXRwdXRbaV0gPSAnPG51bS0nKyhpIC0gOTYpKyc+J1xufVxuXG4vLyBGMS1GMjRcbmZvcihpID0gMTEyOyBpIDwgMTM2OyArK2kpIHtcbiAgb3V0cHV0W2ldID0gJ0YnKyhpLTExMSlcbn1cbiIsIlwidXNlIHN0cmljdFwiXG5cbnZhciBFdmVudEVtaXR0ZXIgPSByZXF1aXJlKFwiZXZlbnRzXCIpLkV2ZW50RW1pdHRlclxuICAsIHV0aWwgICAgICAgICA9IHJlcXVpcmUoXCJ1dGlsXCIpXG4gICwgZG9tcmVhZHkgICAgID0gcmVxdWlyZShcImRvbXJlYWR5XCIpXG4gICwgdmtleSAgICAgICAgID0gcmVxdWlyZShcInZrZXlcIilcbiAgLCBpbnZlcnQgICAgICAgPSByZXF1aXJlKFwiaW52ZXJ0LWhhc2hcIilcbiAgLCB1bmlxICAgICAgICAgPSByZXF1aXJlKFwidW5pcVwiKVxuICAsIGJzZWFyY2ggICAgICA9IHJlcXVpcmUoXCJiaW5hcnktc2VhcmNoLWJvdW5kc1wiKVxuICAsIGlvdGEgICAgICAgICA9IHJlcXVpcmUoXCJpb3RhLWFycmF5XCIpXG4gICwgbWluICAgICAgICAgID0gTWF0aC5taW5cblxuLy9Ccm93c2VyIGNvbXBhdGliaWxpdHkgaGFja3NcbnJlcXVpcmUoXCIuL2xpYi9yYWYtcG9seWZpbGwuanNcIilcbnZhciBhZGRNb3VzZVdoZWVsID0gcmVxdWlyZShcIi4vbGliL21vdXNld2hlZWwtcG9seWZpbGwuanNcIilcbnZhciBocnRpbWUgPSByZXF1aXJlKFwiLi9saWIvaHJ0aW1lLXBvbHlmaWxsLmpzXCIpXG5cbi8vUmVtb3ZlIGFuZ2xlIGJyYWNlcyBhbmQgb3RoZXIgdXNlbGVzcyBjcmFwXG52YXIgZmlsdGVyZWRfdmtleSA9IChmdW5jdGlvbigpIHtcbiAgdmFyIHJlc3VsdCA9IG5ldyBBcnJheSgyNTYpXG4gICAgLCBpLCBqLCBrXG4gIGZvcihpPTA7IGk8MjU2OyArK2kpIHtcbiAgICByZXN1bHRbaV0gPSBcIlVOS1wiXG4gIH1cbiAgZm9yKGkgaW4gdmtleSkge1xuICAgIGsgPSB2a2V5W2ldXG4gICAgaWYoay5jaGFyQXQoMCkgPT09ICc8JyAmJiBrLmNoYXJBdChrLmxlbmd0aC0xKSA9PT0gJz4nKSB7XG4gICAgICBrID0gay5zdWJzdHJpbmcoMSwgay5sZW5ndGgtMSlcbiAgICB9XG4gICAgayA9IGsucmVwbGFjZSgvXFxzL2csIFwiLVwiKVxuICAgIHJlc3VsdFtwYXJzZUludChpKV0gPSBrXG4gIH1cbiAgcmV0dXJuIHJlc3VsdFxufSkoKVxuXG4vL0NvbXB1dGUgbWluaW1hbCBjb21tb24gc2V0IG9mIGtleWJvYXJkIGZ1bmN0aW9uc1xudmFyIGtleU5hbWVzID0gdW5pcShPYmplY3Qua2V5cyhpbnZlcnQoZmlsdGVyZWRfdmtleSkpKVxuXG4vL1RyYW5zbGF0ZXMgYSB2aXJ0dWFsIGtleWNvZGUgdG8gYSBub3JtYWxpemVkIGtleWNvZGVcbmZ1bmN0aW9uIHZpcnR1YWxLZXlDb2RlKGtleSkge1xuICByZXR1cm4gYnNlYXJjaC5lcShrZXlOYW1lcywga2V5KVxufVxuXG4vL01hcHMgYSBwaHlzaWNhbCBrZXljb2RlIHRvIGEgbm9ybWFsaXplZCBrZXljb2RlXG5mdW5jdGlvbiBwaHlzaWNhbEtleUNvZGUoa2V5KSB7XG4gIHJldHVybiB2aXJ0dWFsS2V5Q29kZShmaWx0ZXJlZF92a2V5W2tleV0pXG59XG5cbi8vR2FtZSBzaGVsbFxuZnVuY3Rpb24gR2FtZVNoZWxsKCkge1xuICBFdmVudEVtaXR0ZXIuY2FsbCh0aGlzKVxuICB0aGlzLl9jdXJLZXlTdGF0ZSAgPSBuZXcgQXJyYXkoa2V5TmFtZXMubGVuZ3RoKVxuICB0aGlzLl9wcmVzc0NvdW50ICAgPSBuZXcgQXJyYXkoa2V5TmFtZXMubGVuZ3RoKVxuICB0aGlzLl9yZWxlYXNlQ291bnQgPSBuZXcgQXJyYXkoa2V5TmFtZXMubGVuZ3RoKVxuICBcbiAgdGhpcy5fdGlja0ludGVydmFsID0gbnVsbFxuICB0aGlzLl9yYWZIYW5kbGUgPSBudWxsXG4gIHRoaXMuX3RpY2tSYXRlID0gMFxuICB0aGlzLl9sYXN0VGljayA9IGhydGltZSgpXG4gIHRoaXMuX2ZyYW1lVGltZSA9IDAuMFxuICB0aGlzLl9wYXVzZWQgPSB0cnVlXG4gIHRoaXMuX3dpZHRoID0gMFxuICB0aGlzLl9oZWlnaHQgPSAwXG4gIFxuICB0aGlzLl93YW50RnVsbHNjcmVlbiA9IGZhbHNlXG4gIHRoaXMuX3dhbnRQb2ludGVyTG9jayA9IGZhbHNlXG4gIHRoaXMuX2Z1bGxzY3JlZW5BY3RpdmUgPSBmYWxzZVxuICB0aGlzLl9wb2ludGVyTG9ja0FjdGl2ZSA9IGZhbHNlXG4gIFxuICB0aGlzLl9yZW5kZXIgPSByZW5kZXIuYmluZCh1bmRlZmluZWQsIHRoaXMpXG4gIFxuICBmb3IodmFyIGk9MDsgaTxrZXlOYW1lcy5sZW5ndGg7ICsraSkge1xuICAgIHRoaXMuX2N1cktleVN0YXRlW2ldID0gZmFsc2VcbiAgICB0aGlzLl9wcmVzc0NvdW50W2ldID0gdGhpcy5fcmVsZWFzZUNvdW50W2ldID0gMFxuICB9XG4gIFxuICAvL1B1YmxpYyBtZW1iZXJzXG4gIHRoaXMuZWxlbWVudCA9IG51bGxcbiAgdGhpcy5iaW5kaW5ncyA9IHt9XG4gIHRoaXMuZnJhbWVTa2lwID0gMTAwLjBcbiAgdGhpcy50aWNrQ291bnQgPSAwXG4gIHRoaXMuZnJhbWVDb3VudCA9IDBcbiAgdGhpcy5zdGFydFRpbWUgPSBocnRpbWUoKVxuICB0aGlzLnRpY2tUaW1lID0gdGhpcy5fdGlja1JhdGVcbiAgdGhpcy5mcmFtZVRpbWUgPSAxMC4wXG4gIHRoaXMuc3RpY2t5RnVsbHNjcmVlbiA9IGZhbHNlXG4gIHRoaXMuc3RpY2t5UG9pbnRMb2NrID0gZmFsc2VcbiAgXG4gIC8vU2Nyb2xsIHN0dWZmXG4gIHRoaXMuc2Nyb2xsID0gWzAsMCwwXVxuICAgIFxuICAvL01vdXNlIHN0YXRlXG4gIHRoaXMubW91c2VYID0gMFxuICB0aGlzLm1vdXNlWSA9IDBcbiAgdGhpcy5wcmV2TW91c2VYID0gMFxuICB0aGlzLnByZXZNb3VzZVkgPSAwXG59XG5cbnV0aWwuaW5oZXJpdHMoR2FtZVNoZWxsLCBFdmVudEVtaXR0ZXIpXG5cbnZhciBwcm90byA9IEdhbWVTaGVsbC5wcm90b3R5cGVcblxuLy9CaW5kIGtleW5hbWVzXG5wcm90by5rZXlOYW1lcyA9IGtleU5hbWVzXG5cbi8vQmluZHMgYSB2aXJ0dWFsIGtleWJvYXJkIGV2ZW50IHRvIGEgcGh5c2ljYWwga2V5XG5wcm90by5iaW5kID0gZnVuY3Rpb24odmlydHVhbF9rZXkpIHtcbiAgLy9Mb29rIHVwIHByZXZpb3VzIGtleSBiaW5kaW5nc1xuICB2YXIgYXJyXG4gIGlmKHZpcnR1YWxfa2V5IGluIHRoaXMuYmluZGluZ3MpIHtcbiAgICBhcnIgPSB0aGlzLmJpbmRpbmdzW3ZpcnR1YWxfa2V5XVxuICB9IGVsc2Uge1xuICAgIGFyciA9IFtdXG4gIH1cbiAgLy9BZGQga2V5cyB0byBsaXN0XG4gIHZhciBwaHlzaWNhbF9rZXlcbiAgZm9yKHZhciBpPTEsIG49YXJndW1lbnRzLmxlbmd0aDsgaTxuOyArK2kpIHtcbiAgICBwaHlzaWNhbF9rZXkgPSBhcmd1bWVudHNbaV1cbiAgICBpZih2aXJ0dWFsS2V5Q29kZShwaHlzaWNhbF9rZXkpID49IDApIHtcbiAgICAgIGFyci5wdXNoKHBoeXNpY2FsX2tleSlcbiAgICB9IGVsc2UgaWYocGh5c2ljYWxfa2V5IGluIHRoaXMuYmluZGluZ3MpIHtcbiAgICAgIHZhciBrZXliaW5kcyA9IHRoaXMuYmluZGluZ3NbcGh5c2ljYWxfa2V5XVxuICAgICAgZm9yKHZhciBqPTA7IGo8a2V5YmluZHMubGVuZ3RoOyArK2opIHtcbiAgICAgICAgYXJyLnB1c2goa2V5YmluZHNbal0pXG4gICAgICB9XG4gICAgfVxuICB9XG4gIC8vUmVtb3ZlIGFueSBkdXBsaWNhdGUga2V5c1xuICBhcnIgPSB1bmlxKGFycilcbiAgaWYoYXJyLmxlbmd0aCA+IDApIHtcbiAgICB0aGlzLmJpbmRpbmdzW3ZpcnR1YWxfa2V5XSA9IGFyclxuICB9XG4gIHRoaXMuZW1pdCgnYmluZCcsIHZpcnR1YWxfa2V5LCBhcnIpXG59XG5cbi8vVW5iaW5kcyBhIHZpcnR1YWwga2V5Ym9hcmQgZXZlbnRcbnByb3RvLnVuYmluZCA9IGZ1bmN0aW9uKHZpcnR1YWxfa2V5KSB7XG4gIGlmKHZpcnR1YWxfa2V5IGluIHRoaXMuYmluZGluZ3MpIHtcbiAgICBkZWxldGUgdGhpcy5iaW5kaW5nc1t2aXJ0dWFsX2tleV1cbiAgfVxuICB0aGlzLmVtaXQoJ3VuYmluZCcsIHZpcnR1YWxfa2V5KVxufVxuXG4vL0NoZWNrcyBpZiBhIGtleSBpcyBzZXQgaW4gYSBnaXZlbiBzdGF0ZVxuZnVuY3Rpb24gbG9va3VwS2V5KHN0YXRlLCBiaW5kaW5ncywga2V5KSB7XG4gIGlmKGtleSBpbiBiaW5kaW5ncykge1xuICAgIHZhciBhcnIgPSBiaW5kaW5nc1trZXldXG4gICAgZm9yKHZhciBpPTAsIG49YXJyLmxlbmd0aDsgaTxuOyArK2kpIHtcbiAgICAgIGlmKHN0YXRlW3ZpcnR1YWxLZXlDb2RlKGFycltpXSldKSB7XG4gICAgICAgIHJldHVybiB0cnVlXG4gICAgICB9XG4gICAgfVxuICAgIHJldHVybiBmYWxzZVxuICB9XG4gIHZhciBrYyA9IHZpcnR1YWxLZXlDb2RlKGtleSlcbiAgaWYoa2MgPj0gMCkge1xuICAgIHJldHVybiBzdGF0ZVtrY11cbiAgfVxuICByZXR1cm4gZmFsc2Vcbn1cblxuLy9DaGVja3MgaWYgYSBrZXkgaXMgc2V0IGluIGEgZ2l2ZW4gc3RhdGVcbmZ1bmN0aW9uIGxvb2t1cENvdW50KHN0YXRlLCBiaW5kaW5ncywga2V5KSB7XG4gIGlmKGtleSBpbiBiaW5kaW5ncykge1xuICAgIHZhciBhcnIgPSBiaW5kaW5nc1trZXldLCByID0gMFxuICAgIGZvcih2YXIgaT0wLCBuPWFyci5sZW5ndGg7IGk8bjsgKytpKSB7XG4gICAgICByICs9IHN0YXRlW3ZpcnR1YWxLZXlDb2RlKGFycltpXSldXG4gICAgfVxuICAgIHJldHVybiByXG4gIH1cbiAgdmFyIGtjID0gdmlydHVhbEtleUNvZGUoa2V5KVxuICBpZihrYyA+PSAwKSB7XG4gICAgcmV0dXJuIHN0YXRlW2tjXVxuICB9XG4gIHJldHVybiAwXG59XG5cbi8vQ2hlY2tzIGlmIGEga2V5IChlaXRoZXIgcGh5c2ljYWwgb3IgdmlydHVhbCkgaXMgY3VycmVudGx5IGhlbGQgZG93blxucHJvdG8uZG93biA9IGZ1bmN0aW9uKGtleSkge1xuICByZXR1cm4gbG9va3VwS2V5KHRoaXMuX2N1cktleVN0YXRlLCB0aGlzLmJpbmRpbmdzLCBrZXkpXG59XG5cbi8vQ2hlY2tzIGlmIGEga2V5IHdhcyBldmVyIGRvd25cbnByb3RvLndhc0Rvd24gPSBmdW5jdGlvbihrZXkpIHtcbiAgcmV0dXJuIHRoaXMuZG93bihrZXkpIHx8ICEhdGhpcy5wcmVzcyhrZXkpXG59XG5cbi8vT3Bwb3NpdGUgb2YgZG93blxucHJvdG8udXAgPSBmdW5jdGlvbihrZXkpIHtcbiAgcmV0dXJuICF0aGlzLmRvd24oa2V5KVxufVxuXG4vL0NoZWNrcyBpZiBhIGtleSB3YXMgcmVsZWFzZWQgZHVyaW5nIHByZXZpb3VzIGZyYW1lXG5wcm90by53YXNVcCA9IGZ1bmN0aW9uKGtleSkge1xuICByZXR1cm4gdGhpcy51cChrZXkpIHx8ICEhdGhpcy5yZWxlYXNlKGtleSlcbn1cblxuLy9SZXR1cm5zIHRoZSBudW1iZXIgb2YgdGltZXMgYSBrZXkgd2FzIHByZXNzZWQgc2luY2UgbGFzdCB0aWNrXG5wcm90by5wcmVzcyA9IGZ1bmN0aW9uKGtleSkge1xuICByZXR1cm4gbG9va3VwQ291bnQodGhpcy5fcHJlc3NDb3VudCwgdGhpcy5iaW5kaW5ncywga2V5KVxufVxuXG4vL1JldHVybnMgdGhlIG51bWJlciBvZiB0aW1lcyBhIGtleSB3YXMgcmVsZWFzZWQgc2luY2UgbGFzdCB0aWNrXG5wcm90by5yZWxlYXNlID0gZnVuY3Rpb24oa2V5KSB7XG4gIHJldHVybiBsb29rdXBDb3VudCh0aGlzLl9yZWxlYXNlQ291bnQsIHRoaXMuYmluZGluZ3MsIGtleSlcbn1cblxuLy9QYXVzZS91bnBhdXNlIHRoZSBnYW1lIGxvb3Bcbk9iamVjdC5kZWZpbmVQcm9wZXJ0eShwcm90bywgXCJwYXVzZWRcIiwge1xuICBnZXQ6IGZ1bmN0aW9uKCkge1xuICAgIHJldHVybiB0aGlzLl9wYXVzZWRcbiAgfSxcbiAgc2V0OiBmdW5jdGlvbihzdGF0ZSkge1xuICAgIHZhciBucyA9ICEhc3RhdGVcbiAgICBpZihucyAhPT0gdGhpcy5fcGF1c2VkKSB7XG4gICAgICBpZighdGhpcy5fcGF1c2VkKSB7XG4gICAgICAgIHRoaXMuX3BhdXNlZCA9IHRydWVcbiAgICAgICAgdGhpcy5fZnJhbWVUaW1lID0gbWluKDEuMCwgKGhydGltZSgpIC0gdGhpcy5fbGFzdFRpY2spIC8gdGhpcy5fdGlja1JhdGUpXG4gICAgICAgIGNsZWFySW50ZXJ2YWwodGhpcy5fdGlja0ludGVydmFsKVxuICAgICAgICAvL2NhbmNlbEFuaW1hdGlvbkZyYW1lKHRoaXMuX3JhZkhhbmRsZSlcbiAgICAgIH0gZWxzZSB7XG4gICAgICAgIHRoaXMuX3BhdXNlZCA9IGZhbHNlXG4gICAgICAgIHRoaXMuX2xhc3RUaWNrID0gaHJ0aW1lKCkgLSBNYXRoLmZsb29yKHRoaXMuX2ZyYW1lVGltZSAqIHRoaXMuX3RpY2tSYXRlKVxuICAgICAgICB0aGlzLl90aWNrSW50ZXJ2YWwgPSBzZXRJbnRlcnZhbCh0aWNrLCB0aGlzLl90aWNrUmF0ZSwgdGhpcylcbiAgICAgICAgdGhpcy5fcmFmSGFuZGxlID0gcmVxdWVzdEFuaW1hdGlvbkZyYW1lKHRoaXMuX3JlbmRlcilcbiAgICAgIH1cbiAgICB9XG4gIH1cbn0pXG5cbi8vRnVsbHNjcmVlbiBzdGF0ZSB0b2dnbGVcblxuZnVuY3Rpb24gdHJ5RnVsbHNjcmVlbihzaGVsbCkge1xuICAvL1JlcXVlc3QgZnVsbCBzY3JlZW5cbiAgdmFyIGVsZW0gPSBzaGVsbC5lbGVtZW50XG4gIFxuICBpZihzaGVsbC5fd2FudEZ1bGxzY3JlZW4gJiYgIXNoZWxsLl9mdWxsc2NyZWVuQWN0aXZlKSB7XG4gICAgdmFyIGZzID0gZWxlbS5yZXF1ZXN0RnVsbHNjcmVlbiB8fFxuICAgICAgICAgICAgIGVsZW0ucmVxdWVzdEZ1bGxTY3JlZW4gfHxcbiAgICAgICAgICAgICBlbGVtLndlYmtpdFJlcXVlc3RGdWxsc2NyZWVuIHx8XG4gICAgICAgICAgICAgZWxlbS53ZWJraXRSZXF1ZXN0RnVsbFNjcmVlbiB8fFxuICAgICAgICAgICAgIGVsZW0ubW96UmVxdWVzdEZ1bGxzY3JlZW4gfHxcbiAgICAgICAgICAgICBlbGVtLm1velJlcXVlc3RGdWxsU2NyZWVuIHx8XG4gICAgICAgICAgICAgZnVuY3Rpb24oKSB7fVxuICAgIGZzLmNhbGwoZWxlbSlcbiAgfVxuICBpZihzaGVsbC5fd2FudFBvaW50ZXJMb2NrICYmICFzaGVsbC5fcG9pbnRlckxvY2tBY3RpdmUpIHtcbiAgICB2YXIgcGwgPSAgZWxlbS5yZXF1ZXN0UG9pbnRlckxvY2sgfHxcbiAgICAgICAgICAgICAgZWxlbS53ZWJraXRSZXF1ZXN0UG9pbnRlckxvY2sgfHxcbiAgICAgICAgICAgICAgZWxlbS5tb3pSZXF1ZXN0UG9pbnRlckxvY2sgfHxcbiAgICAgICAgICAgICAgZWxlbS5tc1JlcXVlc3RQb2ludGVyTG9jayB8fFxuICAgICAgICAgICAgICBlbGVtLm9SZXF1ZXN0UG9pbnRlckxvY2sgfHxcbiAgICAgICAgICAgICAgZnVuY3Rpb24oKSB7fVxuICAgIHBsLmNhbGwoZWxlbSlcbiAgfVxufVxuXG52YXIgY2FuY2VsRnVsbHNjcmVlbiA9IGRvY3VtZW50LmV4aXRGdWxsc2NyZWVuIHx8XG4gICAgICAgICAgICAgICAgICAgICAgIGRvY3VtZW50LmNhbmNlbEZ1bGxzY3JlZW4gfHwgIC8vV2h5IGNhbiBubyBvbmUgYWdyZWUgb24gdGhpcz9cbiAgICAgICAgICAgICAgICAgICAgICAgZG9jdW1lbnQuY2FuY2VsRnVsbFNjcmVlbiB8fFxuICAgICAgICAgICAgICAgICAgICAgICBkb2N1bWVudC53ZWJraXRDYW5jZWxGdWxsc2NyZWVuIHx8XG4gICAgICAgICAgICAgICAgICAgICAgIGRvY3VtZW50LndlYmtpdENhbmNlbEZ1bGxTY3JlZW4gfHxcbiAgICAgICAgICAgICAgICAgICAgICAgZG9jdW1lbnQubW96Q2FuY2VsRnVsbHNjcmVlbiB8fFxuICAgICAgICAgICAgICAgICAgICAgICBkb2N1bWVudC5tb3pDYW5jZWxGdWxsU2NyZWVuIHx8XG4gICAgICAgICAgICAgICAgICAgICAgIGZ1bmN0aW9uKCl7fVxuXG5PYmplY3QuZGVmaW5lUHJvcGVydHkocHJvdG8sIFwiZnVsbHNjcmVlblwiLCB7XG4gIGdldDogZnVuY3Rpb24oKSB7XG4gICAgcmV0dXJuIHRoaXMuX2Z1bGxzY3JlZW5BY3RpdmVcbiAgfSxcbiAgc2V0OiBmdW5jdGlvbihzdGF0ZSkge1xuICAgIHZhciBucyA9ICEhc3RhdGVcbiAgICBpZighbnMpIHtcbiAgICAgIHRoaXMuX3dhbnRGdWxsc2NyZWVuID0gZmFsc2VcbiAgICAgIGNhbmNlbEZ1bGxzY3JlZW4uY2FsbChkb2N1bWVudClcbiAgICB9IGVsc2Uge1xuICAgICAgdGhpcy5fd2FudEZ1bGxzY3JlZW4gPSB0cnVlXG4gICAgICB0cnlGdWxsc2NyZWVuKHRoaXMpXG4gICAgfVxuICAgIHJldHVybiB0aGlzLl9mdWxsc2NyZWVuQWN0aXZlXG4gIH1cbn0pXG5cbmZ1bmN0aW9uIGhhbmRsZUZ1bGxzY3JlZW4oc2hlbGwpIHtcbiAgc2hlbGwuX2Z1bGxzY3JlZW5BY3RpdmUgPSBkb2N1bWVudC5mdWxsc2NyZWVuIHx8XG4gICAgICAgICAgICAgICAgICAgICAgICAgICAgZG9jdW1lbnQubW96RnVsbFNjcmVlbiB8fFxuICAgICAgICAgICAgICAgICAgICAgICAgICAgIGRvY3VtZW50LndlYmtpdElzRnVsbFNjcmVlbiB8fFxuICAgICAgICAgICAgICAgICAgICAgICAgICAgIGZhbHNlXG4gIGlmKCFzaGVsbC5zdGlja3lGdWxsc2NyZWVuICYmIHNoZWxsLl9mdWxsc2NyZWVuQWN0aXZlKSB7XG4gICAgc2hlbGwuX3dhbnRGdWxsc2NyZWVuID0gZmFsc2VcbiAgfVxufVxuXG4vL1BvaW50ZXIgbG9jayBzdGF0ZSB0b2dnbGVcbnZhciBleGl0UG9pbnRlckxvY2sgPSBkb2N1bWVudC5leGl0UG9pbnRlckxvY2sgfHxcbiAgICAgICAgICAgICAgICAgICAgICBkb2N1bWVudC53ZWJraXRFeGl0UG9pbnRlckxvY2sgfHxcbiAgICAgICAgICAgICAgICAgICAgICBkb2N1bWVudC5tb3pFeGl0UG9pbnRlckxvY2sgfHxcbiAgICAgICAgICAgICAgICAgICAgICBmdW5jdGlvbigpIHt9XG5cbk9iamVjdC5kZWZpbmVQcm9wZXJ0eShwcm90bywgXCJwb2ludGVyTG9ja1wiLCB7XG4gIGdldDogZnVuY3Rpb24oKSB7XG4gICAgcmV0dXJuIHRoaXMuX3BvaW50ZXJMb2NrQWN0aXZlXG4gIH0sXG4gIHNldDogZnVuY3Rpb24oc3RhdGUpIHtcbiAgICB2YXIgbnMgPSAhIXN0YXRlXG4gICAgaWYoIW5zKSB7XG4gICAgICB0aGlzLl93YW50UG9pbnRlckxvY2sgPSBmYWxzZVxuICAgICAgZXhpdFBvaW50ZXJMb2NrLmNhbGwoZG9jdW1lbnQpXG4gICAgfSBlbHNlIHtcbiAgICAgIHRoaXMuX3dhbnRQb2ludGVyTG9jayA9IHRydWVcbiAgICAgIHRyeUZ1bGxzY3JlZW4odGhpcylcbiAgICB9XG4gICAgcmV0dXJuIHRoaXMuX3BvaW50ZXJMb2NrQWN0aXZlXG4gIH1cbn0pXG5cbmZ1bmN0aW9uIGhhbmRsZVBvaW50ZXJMb2NrQ2hhbmdlKHNoZWxsLCBldmVudCkge1xuICBzaGVsbC5fcG9pbnRlckxvY2tBY3RpdmUgPSBzaGVsbC5lbGVtZW50ID09PSAoXG4gICAgICBkb2N1bWVudC5wb2ludGVyTG9ja0VsZW1lbnQgfHxcbiAgICAgIGRvY3VtZW50Lm1velBvaW50ZXJMb2NrRWxlbWVudCB8fFxuICAgICAgZG9jdW1lbnQud2Via2l0UG9pbnRlckxvY2tFbGVtZW50IHx8XG4gICAgICBudWxsKVxuICBpZighc2hlbGwuc3RpY2t5UG9pbnRlckxvY2sgJiYgc2hlbGwuX3BvaW50ZXJMb2NrQWN0aXZlKSB7XG4gICAgc2hlbGwuX3dhbnRQb2ludGVyTG9jayA9IGZhbHNlXG4gIH1cbn1cblxuLy9XaWR0aCBhbmQgaGVpZ2h0XG5PYmplY3QuZGVmaW5lUHJvcGVydHkocHJvdG8sIFwid2lkdGhcIiwge1xuICBnZXQ6IGZ1bmN0aW9uKCkge1xuICAgIHJldHVybiB0aGlzLmVsZW1lbnQuY2xpZW50V2lkdGhcbiAgfVxufSlcbk9iamVjdC5kZWZpbmVQcm9wZXJ0eShwcm90bywgXCJoZWlnaHRcIiwge1xuICBnZXQ6IGZ1bmN0aW9uKCkge1xuICAgIHJldHVybiB0aGlzLmVsZW1lbnQuY2xpZW50SGVpZ2h0XG4gIH1cbn0pXG5cbi8vU2V0IGtleSBzdGF0ZVxuZnVuY3Rpb24gc2V0S2V5U3RhdGUoc2hlbGwsIGtleSwgc3RhdGUpIHtcbiAgdmFyIHBzID0gc2hlbGwuX2N1cktleVN0YXRlW2tleV1cbiAgaWYocHMgIT09IHN0YXRlKSB7XG4gICAgaWYoc3RhdGUpIHtcbiAgICAgIHNoZWxsLl9wcmVzc0NvdW50W2tleV0rK1xuICAgIH0gZWxzZSB7XG4gICAgICBzaGVsbC5fcmVsZWFzZUNvdW50W2tleV0rK1xuICAgIH1cbiAgICBzaGVsbC5fY3VyS2V5U3RhdGVba2V5XSA9IHN0YXRlXG4gIH1cbn1cblxuLy9UaWNrcyB0aGUgZ2FtZSBzdGF0ZSBvbmUgdXBkYXRlXG5mdW5jdGlvbiB0aWNrKHNoZWxsKSB7XG4gIHZhciBza2lwID0gaHJ0aW1lKCkgKyBzaGVsbC5mcmFtZVNraXBcbiAgICAsIHBDb3VudCA9IHNoZWxsLl9wcmVzc0NvdW50XG4gICAgLCByQ291bnQgPSBzaGVsbC5fcmVsZWFzZUNvdW50XG4gICAgLCBpLCBzLCB0XG4gICAgLCB0ciA9IHNoZWxsLl90aWNrUmF0ZVxuICAgICwgbiA9IGtleU5hbWVzLmxlbmd0aFxuICB3aGlsZSghc2hlbGwuX3BhdXNlZCAmJlxuICAgICAgICBocnRpbWUoKSA+PSBzaGVsbC5fbGFzdFRpY2sgKyB0cikge1xuICAgIFxuICAgIC8vU2tpcCBmcmFtZXMgaWYgd2UgYXJlIG92ZXIgYnVkZ2V0XG4gICAgaWYoaHJ0aW1lKCkgPiBza2lwKSB7XG4gICAgICBzaGVsbC5fbGFzdFRpY2sgPSBocnRpbWUoKSArIHRyXG4gICAgICByZXR1cm5cbiAgICB9XG4gICAgXG4gICAgLy9UaWNrIHRoZSBnYW1lXG4gICAgcyA9IGhydGltZSgpXG4gICAgc2hlbGwuZW1pdChcInRpY2tcIilcbiAgICB0ID0gaHJ0aW1lKClcbiAgICBzaGVsbC50aWNrVGltZSA9IHQgLSBzXG4gICAgXG4gICAgLy9VcGRhdGUgY291bnRlcnMgYW5kIHRpbWVcbiAgICArK3NoZWxsLnRpY2tDb3VudFxuICAgIHNoZWxsLl9sYXN0VGljayArPSB0clxuICAgIFxuICAgIC8vU2hpZnQgaW5wdXQgc3RhdGVcbiAgICBmb3IoaT0wOyBpPG47ICsraSkge1xuICAgICAgcENvdW50W2ldID0gckNvdW50W2ldID0gMFxuICAgIH1cbiAgICBpZihzaGVsbC5fcG9pbnRlckxvY2tBY3RpdmUpIHtcbiAgICAgIHNoZWxsLnByZXZNb3VzZVggPSBzaGVsbC5tb3VzZVggPSBzaGVsbC53aWR0aD4+MVxuICAgICAgc2hlbGwucHJldk1vdXNlWSA9IHNoZWxsLm1vdXNlWSA9IHNoZWxsLmhlaWdodD4+MVxuICAgIH0gZWxzZSB7XG4gICAgICBzaGVsbC5wcmV2TW91c2VYID0gc2hlbGwubW91c2VYXG4gICAgICBzaGVsbC5wcmV2TW91c2VZID0gc2hlbGwubW91c2VZXG4gICAgfVxuICAgIHNoZWxsLnNjcm9sbFswXSA9IHNoZWxsLnNjcm9sbFsxXSA9IHNoZWxsLnNjcm9sbFsyXSA9IDBcbiAgfVxufVxuXG4vL1JlbmRlciBzdHVmZlxuZnVuY3Rpb24gcmVuZGVyKHNoZWxsKSB7XG5cbiAgLy9SZXF1ZXN0IG5leHQgZnJhbWVcbiAgc2hlbGwuX3JhZkhhbmRsZSA9IHJlcXVlc3RBbmltYXRpb25GcmFtZShzaGVsbC5fcmVuZGVyKVxuXG4gIC8vVGljayB0aGUgc2hlbGxcbiAgdGljayhzaGVsbClcbiAgXG4gIC8vQ29tcHV0ZSBmcmFtZSB0aW1lXG4gIHZhciBkdFxuICBpZihzaGVsbC5fcGF1c2VkKSB7XG4gICAgZHQgPSBzaGVsbC5fZnJhbWVUaW1lXG4gIH0gZWxzZSB7XG4gICAgZHQgPSBtaW4oMS4wLCAoaHJ0aW1lKCkgLSBzaGVsbC5fbGFzdFRpY2spIC8gc2hlbGwuX3RpY2tSYXRlKVxuICB9XG4gIFxuICAvL0RyYXcgYSBmcmFtZVxuICArK3NoZWxsLmZyYW1lQ291bnRcbiAgdmFyIHMgPSBocnRpbWUoKVxuICBzaGVsbC5lbWl0KFwicmVuZGVyXCIsIGR0KVxuICB2YXIgdCA9IGhydGltZSgpXG4gIHNoZWxsLmZyYW1lVGltZSA9IHQgLSBzXG4gIFxufVxuXG5mdW5jdGlvbiBpc0ZvY3VzZWQoc2hlbGwpIHtcbiAgcmV0dXJuIChkb2N1bWVudC5hY3RpdmVFbGVtZW50ID09PSBkb2N1bWVudC5ib2R5KSB8fFxuICAgICAgICAgKGRvY3VtZW50LmFjdGl2ZUVsZW1lbnQgPT09IHNoZWxsLmVsZW1lbnQpXG59XG5cbi8vU2V0IGtleSB1cFxuZnVuY3Rpb24gaGFuZGxlS2V5VXAoc2hlbGwsIGV2KSB7XG4gIGV2LnByZXZlbnREZWZhdWx0KClcbiAgdmFyIGtjID0gcGh5c2ljYWxLZXlDb2RlKGV2LmtleUNvZGUgfHwgZXYuY2hhciB8fCBldi53aGljaCB8fCBldi5jaGFyQ29kZSlcbiAgaWYoa2MgPj0gMCkge1xuICAgIHNldEtleVN0YXRlKHNoZWxsLCBrYywgZmFsc2UpXG4gIH1cbn1cblxuLy9TZXQga2V5IGRvd25cbmZ1bmN0aW9uIGhhbmRsZUtleURvd24oc2hlbGwsIGV2KSB7XG4gIGlmKCFpc0ZvY3VzZWQoc2hlbGwpKSB7XG4gICAgcmV0dXJuXG4gIH1cbiAgaWYoZXYubWV0YUtleSkge1xuICAgIC8vSGFjazogQ2xlYXIga2V5IHN0YXRlIHdoZW4gbWV0YSBnZXRzIHByZXNzZWQgdG8gcHJldmVudCBrZXlzIHN0aWNraW5nXG4gICAgaGFuZGxlQmx1cihzaGVsbCwgZXYpXG4gIH0gZWxzZSB7XG4gICAgZXYucHJldmVudERlZmF1bHQoKVxuICAgIHZhciBrYyA9IHBoeXNpY2FsS2V5Q29kZShldi5rZXlDb2RlIHx8IGV2LmNoYXIgfHwgZXYud2hpY2ggfHwgZXYuY2hhckNvZGUpXG4gICAgaWYoa2MgPj0gMCkge1xuICAgICAgc2V0S2V5U3RhdGUoc2hlbGwsIGtjLCB0cnVlKVxuICAgIH1cbiAgfVxufVxuXG4vL01vdXNlIGV2ZW50cyBhcmUgcmVhbGx5IGFubm95aW5nXG52YXIgbW91c2VDb2RlcyA9IGlvdGEoMzIpLm1hcChmdW5jdGlvbihuKSB7XG4gIHJldHVybiB2aXJ0dWFsS2V5Q29kZShcIm1vdXNlLVwiICsgKG4rMSkpXG59KVxuXG5mdW5jdGlvbiBzZXRNb3VzZUJ1dHRvbnMoc2hlbGwsIGJ1dHRvbnMpIHtcbiAgZm9yKHZhciBpPTA7IGk8MzI7ICsraSkge1xuICAgIHNldEtleVN0YXRlKHNoZWxsLCBtb3VzZUNvZGVzW2ldLCAhIShidXR0b25zICYgKDE8PGkpKSlcbiAgfVxufVxuXG5mdW5jdGlvbiBoYW5kbGVNb3VzZU1vdmUoc2hlbGwsIGV2KSB7XG4gIGlmKHNoZWxsLl9wb2ludGVyTG9ja0FjdGl2ZSkge1xuICAgIHZhciBtb3ZlbWVudFggPSBldi5tb3ZlbWVudFggICAgICAgfHxcbiAgICAgICAgICAgICAgICAgICAgZXYubW96TW92ZW1lbnRYICAgIHx8XG4gICAgICAgICAgICAgICAgICAgIGV2LndlYmtpdE1vdmVtZW50WCB8fFxuICAgICAgICAgICAgICAgICAgICAwLFxuICAgICAgICBtb3ZlbWVudFkgPSBldi5tb3ZlbWVudFkgICAgICAgfHxcbiAgICAgICAgICAgICAgICAgICAgZXYubW96TW92ZW1lbnRZICAgIHx8XG4gICAgICAgICAgICAgICAgICAgIGV2LndlYmtpdE1vdmVtZW50WSB8fFxuICAgICAgICAgICAgICAgICAgICAwXG4gICAgc2hlbGwubW91c2VYICs9IG1vdmVtZW50WFxuICAgIHNoZWxsLm1vdXNlWSArPSBtb3ZlbWVudFlcbiAgfSBlbHNlIHtcbiAgICBzaGVsbC5tb3VzZVggPSBldi5jbGllbnRYIC0gc2hlbGwuZWxlbWVudC5vZmZzZXRMZWZ0XG4gICAgc2hlbGwubW91c2VZID0gZXYuY2xpZW50WSAtIHNoZWxsLmVsZW1lbnQub2Zmc2V0VG9wXG4gIH1cbiAgcmV0dXJuIGZhbHNlXG59XG5cbmZ1bmN0aW9uIGhhbmRsZU1vdXNlRG93bihzaGVsbCwgZXYpIHtcbiAgc2V0S2V5U3RhdGUoc2hlbGwsIG1vdXNlQ29kZXNbZXYuYnV0dG9uXSwgdHJ1ZSlcbiAgcmV0dXJuIGZhbHNlXG59XG5cbmZ1bmN0aW9uIGhhbmRsZU1vdXNlVXAoc2hlbGwsIGV2KSB7XG4gIHNldEtleVN0YXRlKHNoZWxsLCBtb3VzZUNvZGVzW2V2LmJ1dHRvbl0sIGZhbHNlKVxuICByZXR1cm4gZmFsc2Vcbn1cblxuZnVuY3Rpb24gaGFuZGxlTW91c2VFbnRlcihzaGVsbCwgZXYpIHtcbiAgaWYoc2hlbGwuX3BvaW50ZXJMb2NrQWN0aXZlKSB7XG4gICAgc2hlbGwucHJldk1vdXNlWCA9IHNoZWxsLm1vdXNlWCA9IHNoZWxsLndpZHRoPj4xXG4gICAgc2hlbGwucHJldk1vdXNlWSA9IHNoZWxsLm1vdXNlWSA9IHNoZWxsLmhlaWdodD4+MVxuICB9IGVsc2Uge1xuICAgIHNoZWxsLnByZXZNb3VzZVggPSBzaGVsbC5tb3VzZVggPSBldi5jbGllbnRYIC0gc2hlbGwuZWxlbWVudC5vZmZzZXRMZWZ0XG4gICAgc2hlbGwucHJldk1vdXNlWSA9IHNoZWxsLm1vdXNlWSA9IGV2LmNsaWVudFkgLSBzaGVsbC5lbGVtZW50Lm9mZnNldFRvcFxuICB9XG4gIHJldHVybiBmYWxzZVxufVxuXG5mdW5jdGlvbiBoYW5kbGVNb3VzZUxlYXZlKHNoZWxsLCBldikge1xuICBzZXRNb3VzZUJ1dHRvbnMoc2hlbGwsIDApXG4gIHJldHVybiBmYWxzZVxufVxuXG4vL0hhbmRsZSBtb3VzZSB3aGVlbCBldmVudHNcbmZ1bmN0aW9uIGhhbmRsZU1vdXNlV2hlZWwoc2hlbGwsIGV2KSB7XG4gIHZhciBzY2FsZSA9IDFcbiAgc3dpdGNoKGV2LmRlbHRhTW9kZSkge1xuICAgIGNhc2UgMDogLy9QaXhlbFxuICAgICAgc2NhbGUgPSAxXG4gICAgYnJlYWtcbiAgICBjYXNlIDE6IC8vTGluZVxuICAgICAgc2NhbGUgPSAxMlxuICAgIGJyZWFrXG4gICAgY2FzZSAyOiAvL1BhZ2VcbiAgICAgICBzY2FsZSA9IHNoZWxsLmhlaWdodFxuICAgIGJyZWFrXG4gIH1cbiAgLy9BZGQgc2Nyb2xsXG4gIHNoZWxsLnNjcm9sbFswXSArPSAgZXYuZGVsdGFYICogc2NhbGVcbiAgc2hlbGwuc2Nyb2xsWzFdICs9ICBldi5kZWx0YVkgKiBzY2FsZVxuICBzaGVsbC5zY3JvbGxbMl0gKz0gKGV2LmRlbHRhWiAqIHNjYWxlKXx8MC4wXG4gIHJldHVybiBmYWxzZVxufVxuXG5mdW5jdGlvbiBoYW5kbGVDb250ZXhNZW51KHNoZWxsLCBldikge1xuICByZXR1cm4gZmFsc2Vcbn1cblxuZnVuY3Rpb24gaGFuZGxlQmx1cihzaGVsbCwgZXYpIHtcbiAgdmFyIG4gPSBrZXlOYW1lcy5sZW5ndGhcbiAgICAsIGMgPSBzaGVsbC5fY3VyS2V5U3RhdGVcbiAgICAsIHIgPSBzaGVsbC5fcmVsZWFzZUNvdW50XG4gICAgLCBpXG4gIGZvcihpPTA7IGk8bjsgKytpKSB7XG4gICAgaWYoY1tpXSkge1xuICAgICAgKytyW2ldXG4gICAgfVxuICAgIGNbaV0gPSBmYWxzZVxuICB9XG4gIHJldHVybiBmYWxzZVxufVxuXG5mdW5jdGlvbiBoYW5kbGVSZXNpemVFbGVtZW50KHNoZWxsLCBldikge1xuICB2YXIgdyA9IHNoZWxsLmVsZW1lbnQuY2xpZW50V2lkdGh8MFxuICB2YXIgaCA9IHNoZWxsLmVsZW1lbnQuY2xpZW50SGVpZ2h0fDBcbiAgaWYoKHcgIT09IHNoZWxsLl93aWR0aCkgfHwgKGggIT09IHNoZWxsLl9oZWlnaHQpKSB7XG4gICAgc2hlbGwuX3dpZHRoID0gd1xuICAgIHNoZWxsLl9oZWlnaHQgPSBoXG4gICAgc2hlbGwuZW1pdChcInJlc2l6ZVwiLCB3LCBoKVxuICB9XG59XG5cbmZ1bmN0aW9uIG1ha2VEZWZhdWx0Q29udGFpbmVyKCkge1xuICB2YXIgY29udGFpbmVyID0gZG9jdW1lbnQuY3JlYXRlRWxlbWVudChcImRpdlwiKVxuICBjb250YWluZXIudGFiaW5kZXggPSAxXG4gIGNvbnRhaW5lci5zdHlsZS5wb3NpdGlvbiA9IFwiYWJzb2x1dGVcIlxuICBjb250YWluZXIuc3R5bGUubGVmdCA9IFwiMHB4XCJcbiAgY29udGFpbmVyLnN0eWxlLnJpZ2h0ID0gXCIwcHhcIlxuICBjb250YWluZXIuc3R5bGUudG9wID0gXCIwcHhcIlxuICBjb250YWluZXIuc3R5bGUuYm90dG9tID0gXCIwcHhcIlxuICBjb250YWluZXIuc3R5bGUuaGVpZ2h0ID0gXCIxMDAlXCJcbiAgY29udGFpbmVyLnN0eWxlLm92ZXJmbG93ID0gXCJoaWRkZW5cIlxuICBkb2N1bWVudC5ib2R5LmFwcGVuZENoaWxkKGNvbnRhaW5lcilcbiAgZG9jdW1lbnQuYm9keS5zdHlsZS5vdmVyZmxvdyA9IFwiaGlkZGVuXCIgLy9QcmV2ZW50IGJvdW5jZVxuICBkb2N1bWVudC5ib2R5LnN0eWxlLmhlaWdodCA9IFwiMTAwJVwiXG4gIHJldHVybiBjb250YWluZXJcbn1cblxuZnVuY3Rpb24gY3JlYXRlU2hlbGwob3B0aW9ucykge1xuICBvcHRpb25zID0gb3B0aW9ucyB8fCB7fVxuICBcbiAgLy9DaGVjayBmdWxsc2NyZWVuIGFuZCBwb2ludGVyIGxvY2sgZmxhZ3NcbiAgdmFyIHVzZUZ1bGxzY3JlZW4gPSAhIW9wdGlvbnMuZnVsbHNjcmVlblxuICB2YXIgdXNlUG9pbnRlckxvY2sgPSB1c2VGdWxsc2NyZWVuXG4gIGlmKHR5cGVvZiBvcHRpb25zLnBvaW50ZXJMb2NrICE9PSB1bmRlZmluZWQpIHtcbiAgICB1c2VQb2ludGVyTG9jayA9ICEhb3B0aW9ucy5wb2ludGVyTG9ja1xuICB9XG4gIFxuICAvL0NyZWF0ZSBpbml0aWFsIHNoZWxsXG4gIHZhciBzaGVsbCA9IG5ldyBHYW1lU2hlbGwoKVxuICBzaGVsbC5fdGlja1JhdGUgPSBvcHRpb25zLnRpY2tSYXRlIHx8IDMwXG4gIHNoZWxsLmZyYW1lU2tpcCA9IG9wdGlvbnMuZnJhbWVTa2lwIHx8IChzaGVsbC5fdGlja1JhdGUrNSkgKiA1XG4gIHNoZWxsLnN0aWNreUZ1bGxzY3JlZW4gPSAhIW9wdGlvbnMuc3RpY2t5RnVsbHNjcmVlbiB8fCAhIW9wdGlvbnMuc3RpY2t5XG4gIHNoZWxsLnN0aWNreVBvaW50ZXJMb2NrID0gISFvcHRpb25zLnN0aWNrUG9pbnRlckxvY2sgfHwgIW9wdGlvbnMuc3RpY2t5XG4gIFxuICAvL1NldCBiaW5kaW5nc1xuICBpZihvcHRpb25zLmJpbmRpbmdzKSB7XG4gICAgc2hlbGwuYmluZGluZ3MgPSBiaW5kaW5nc1xuICB9XG4gIFxuICAvL1dhaXQgZm9yIGRvbSB0byBpbnRpYWlsaXplXG4gIHNldFRpbWVvdXQoZnVuY3Rpb24oKSB7IGRvbXJlYWR5KGZ1bmN0aW9uIGluaXRHYW1lU2hlbGwoKSB7XG4gICAgXG4gICAgLy9SZXRyaWV2ZSBlbGVtZW50XG4gICAgdmFyIGVsZW1lbnQgPSBvcHRpb25zLmVsZW1lbnRcbiAgICBpZih0eXBlb2YgZWxlbWVudCA9PT0gXCJzdHJpbmdcIikge1xuICAgICAgdmFyIGUgPSBkb2N1bWVudC5xdWVyeVNlbGVjdG9yKGVsZW1lbnQpXG4gICAgICBpZighZSkge1xuICAgICAgICBlID0gZG9jdW1lbnQuZ2V0RWxlbWVudEJ5SWQoZWxlbWVudClcbiAgICAgIH1cbiAgICAgIGlmKCFlKSB7XG4gICAgICAgIGUgPSBkb2N1bWVudC5nZXRFbGVtZW50QnlDbGFzcyhlbGVtZW50KVswXVxuICAgICAgfVxuICAgICAgaWYoIWUpIHtcbiAgICAgICAgZSA9IG1ha2VEZWZhdWx0Q29udGFpbmVyKClcbiAgICAgIH1cbiAgICAgIHNoZWxsLmVsZW1lbnQgPSBlXG4gICAgfSBlbHNlIGlmKHR5cGVvZiBlbGVtZW50ID09PSBcIm9iamVjdFwiICYmICEhZWxlbWVudCkge1xuICAgICAgc2hlbGwuZWxlbWVudCA9IGVsZW1lbnRcbiAgICB9IGVsc2UgaWYodHlwZW9mIGVsZW1lbnQgPT09IFwiZnVuY3Rpb25cIikge1xuICAgICAgc2hlbGwuZWxlbWVudCA9IGVsZW1lbnQoKVxuICAgIH0gZWxzZSB7XG4gICAgICBzaGVsbC5lbGVtZW50ID0gbWFrZURlZmF1bHRDb250YWluZXIoKVxuICAgIH1cbiAgICBcbiAgICAvL0Rpc2FibGUgdXNlci1zZWxlY3RcbiAgICBpZihzaGVsbC5lbGVtZW50LnN0eWxlKSB7XG4gICAgICBzaGVsbC5lbGVtZW50LnN0eWxlW1wiLXdlYmtpdC10b3VjaC1jYWxsb3V0XCJdID0gXCJub25lXCJcbiAgICAgIHNoZWxsLmVsZW1lbnQuc3R5bGVbXCItd2Via2l0LXVzZXItc2VsZWN0XCJdID0gXCJub25lXCJcbiAgICAgIHNoZWxsLmVsZW1lbnQuc3R5bGVbXCIta2h0bWwtdXNlci1zZWxlY3RcIl0gPSBcIm5vbmVcIlxuICAgICAgc2hlbGwuZWxlbWVudC5zdHlsZVtcIi1tb3otdXNlci1zZWxlY3RcIl0gPSBcIm5vbmVcIlxuICAgICAgc2hlbGwuZWxlbWVudC5zdHlsZVtcIi1tcy11c2VyLXNlbGVjdFwiXSA9IFwibm9uZVwiXG4gICAgICBzaGVsbC5lbGVtZW50LnN0eWxlW1widXNlci1zZWxlY3RcIl0gPSBcIm5vbmVcIlxuICAgIH1cbiAgICBcbiAgICAvL0hvb2sgcmVzaXplIGhhbmRsZXJcbiAgICBzaGVsbC5fd2lkdGggPSBzaGVsbC5lbGVtZW50LmNsaWVudFdpZHRoXG4gICAgc2hlbGwuX2hlaWdodCA9IHNoZWxsLmVsZW1lbnQuY2xpZW50SGVpZ2h0XG4gICAgdmFyIGhhbmRsZVJlc2l6ZSA9IGhhbmRsZVJlc2l6ZUVsZW1lbnQuYmluZCh1bmRlZmluZWQsIHNoZWxsKVxuICAgIGlmKHR5cGVvZiBNdXRhdGlvbk9ic2VydmVyICE9PSBcInVuZGVmaW5lZFwiKSB7XG4gICAgICB2YXIgb2JzZXJ2ZXIgPSBuZXcgTXV0YXRpb25PYnNlcnZlcihoYW5kbGVSZXNpemUpXG4gICAgICBvYnNlcnZlci5vYnNlcnZlKHNoZWxsLmVsZW1lbnQsIHtcbiAgICAgICAgYXR0cmlidXRlczogdHJ1ZSxcbiAgICAgICAgc3VidHJlZTogdHJ1ZVxuICAgICAgfSlcbiAgICB9IGVsc2Uge1xuICAgICAgc2hlbGwuZWxlbWVudC5hZGRFdmVudExpc3RlbmVyKFwiRE9NU3VidHJlZU1vZGlmaWVkXCIsIGhhbmRsZVJlc2l6ZSwgZmFsc2UpXG4gICAgfVxuICAgIHdpbmRvdy5hZGRFdmVudExpc3RlbmVyKFwicmVzaXplXCIsIGhhbmRsZVJlc2l6ZSwgZmFsc2UpXG4gICAgXG4gICAgLy9Ib29rIGtleWJvYXJkIGxpc3RlbmVyXG4gICAgd2luZG93LmFkZEV2ZW50TGlzdGVuZXIoXCJrZXlkb3duXCIsIGhhbmRsZUtleURvd24uYmluZCh1bmRlZmluZWQsIHNoZWxsKSwgZmFsc2UpXG4gICAgd2luZG93LmFkZEV2ZW50TGlzdGVuZXIoXCJrZXl1cFwiLCBoYW5kbGVLZXlVcC5iaW5kKHVuZGVmaW5lZCwgc2hlbGwpLCBmYWxzZSlcbiAgICBcbiAgICAvL0Rpc2FibGUgcmlnaHQgY2xpY2tcbiAgICBzaGVsbC5lbGVtZW50Lm9uY29udGV4dG1lbnUgPSBoYW5kbGVDb250ZXhNZW51LmJpbmQodW5kZWZpbmVkLCBzaGVsbClcbiAgICBcbiAgICAvL0hvb2sgbW91c2UgbGlzdGVuZXJzXG4gICAgc2hlbGwuZWxlbWVudC5hZGRFdmVudExpc3RlbmVyKFwibW91c2Vkb3duXCIsIGhhbmRsZU1vdXNlRG93bi5iaW5kKHVuZGVmaW5lZCwgc2hlbGwpLCBmYWxzZSlcbiAgICBzaGVsbC5lbGVtZW50LmFkZEV2ZW50TGlzdGVuZXIoXCJtb3VzZXVwXCIsIGhhbmRsZU1vdXNlVXAuYmluZCh1bmRlZmluZWQsIHNoZWxsKSwgZmFsc2UpXG4gICAgc2hlbGwuZWxlbWVudC5hZGRFdmVudExpc3RlbmVyKFwibW91c2Vtb3ZlXCIsIGhhbmRsZU1vdXNlTW92ZS5iaW5kKHVuZGVmaW5lZCwgc2hlbGwpLCBmYWxzZSlcbiAgICBzaGVsbC5lbGVtZW50LmFkZEV2ZW50TGlzdGVuZXIoXCJtb3VzZWVudGVyXCIsIGhhbmRsZU1vdXNlRW50ZXIuYmluZCh1bmRlZmluZWQsIHNoZWxsKSwgZmFsc2UpXG4gICAgXG4gICAgLy9Nb3VzZSBsZWF2ZVxuICAgIHZhciBsZWF2ZSA9IGhhbmRsZU1vdXNlTGVhdmUuYmluZCh1bmRlZmluZWQsIHNoZWxsKVxuICAgIHNoZWxsLmVsZW1lbnQuYWRkRXZlbnRMaXN0ZW5lcihcIm1vdXNlbGVhdmVcIiwgbGVhdmUsIGZhbHNlKVxuICAgIHNoZWxsLmVsZW1lbnQuYWRkRXZlbnRMaXN0ZW5lcihcIm1vdXNlb3V0XCIsIGxlYXZlLCBmYWxzZSlcbiAgICB3aW5kb3cuYWRkRXZlbnRMaXN0ZW5lcihcIm1vdXNlbGVhdmVcIiwgbGVhdmUsIGZhbHNlKVxuICAgIHdpbmRvdy5hZGRFdmVudExpc3RlbmVyKFwibW91c2VvdXRcIiwgbGVhdmUsIGZhbHNlKVxuICAgIFxuICAgIC8vQmx1ciBldmVudCBcbiAgICB2YXIgYmx1ciA9IGhhbmRsZUJsdXIuYmluZCh1bmRlZmluZWQsIHNoZWxsKVxuICAgIHNoZWxsLmVsZW1lbnQuYWRkRXZlbnRMaXN0ZW5lcihcImJsdXJcIiwgYmx1ciwgZmFsc2UpXG4gICAgc2hlbGwuZWxlbWVudC5hZGRFdmVudExpc3RlbmVyKFwiZm9jdXNvdXRcIiwgYmx1ciwgZmFsc2UpXG4gICAgc2hlbGwuZWxlbWVudC5hZGRFdmVudExpc3RlbmVyKFwiZm9jdXNcIiwgYmx1ciwgZmFsc2UpXG4gICAgd2luZG93LmFkZEV2ZW50TGlzdGVuZXIoXCJibHVyXCIsIGJsdXIsIGZhbHNlKVxuICAgIHdpbmRvdy5hZGRFdmVudExpc3RlbmVyKFwiZm9jdXNvdXRcIiwgYmx1ciwgZmFsc2UpXG4gICAgd2luZG93LmFkZEV2ZW50TGlzdGVuZXIoXCJmb2N1c1wiLCBibHVyLCBmYWxzZSlcblxuICAgIC8vTW91c2Ugd2hlZWwgaGFuZGxlclxuICAgIGFkZE1vdXNlV2hlZWwoc2hlbGwuZWxlbWVudCwgaGFuZGxlTW91c2VXaGVlbC5iaW5kKHVuZGVmaW5lZCwgc2hlbGwpLCBmYWxzZSlcblxuICAgIC8vRnVsbHNjcmVlbiBoYW5kbGVyXG4gICAgdmFyIGZ1bGxzY3JlZW5DaGFuZ2UgPSBoYW5kbGVGdWxsc2NyZWVuLmJpbmQodW5kZWZpbmVkLCBzaGVsbClcbiAgICBkb2N1bWVudC5hZGRFdmVudExpc3RlbmVyKFwiZnVsbHNjcmVlbmNoYW5nZVwiLCBmdWxsc2NyZWVuQ2hhbmdlLCBmYWxzZSlcbiAgICBkb2N1bWVudC5hZGRFdmVudExpc3RlbmVyKFwibW96ZnVsbHNjcmVlbmNoYW5nZVwiLCBmdWxsc2NyZWVuQ2hhbmdlLCBmYWxzZSlcbiAgICBkb2N1bWVudC5hZGRFdmVudExpc3RlbmVyKFwid2Via2l0ZnVsbHNjcmVlbmNoYW5nZVwiLCBmdWxsc2NyZWVuQ2hhbmdlLCBmYWxzZSlcblxuICAgIC8vU3R1cGlkIGZ1bGxzY3JlZW4gaGFja1xuICAgIHNoZWxsLmVsZW1lbnQuYWRkRXZlbnRMaXN0ZW5lcihcImNsaWNrXCIsIHRyeUZ1bGxzY3JlZW4uYmluZCh1bmRlZmluZWQsIHNoZWxsKSwgZmFsc2UpXG5cbiAgICAvL1BvaW50ZXIgbG9jayBjaGFuZ2UgaGFuZGxlclxuICAgIHZhciBwb2ludGVyTG9ja0NoYW5nZSA9IGhhbmRsZVBvaW50ZXJMb2NrQ2hhbmdlLmJpbmQodW5kZWZpbmVkLCBzaGVsbClcbiAgICBkb2N1bWVudC5hZGRFdmVudExpc3RlbmVyKFwicG9pbnRlcmxvY2tjaGFuZ2VcIiwgcG9pbnRlckxvY2tDaGFuZ2UsIGZhbHNlKVxuICAgIGRvY3VtZW50LmFkZEV2ZW50TGlzdGVuZXIoXCJtb3pwb2ludGVybG9ja2NoYW5nZVwiLCBwb2ludGVyTG9ja0NoYW5nZSwgZmFsc2UpXG4gICAgZG9jdW1lbnQuYWRkRXZlbnRMaXN0ZW5lcihcIndlYmtpdHBvaW50ZXJsb2NrY2hhbmdlXCIsIHBvaW50ZXJMb2NrQ2hhbmdlLCBmYWxzZSlcbiAgICBkb2N1bWVudC5hZGRFdmVudExpc3RlbmVyKFwicG9pbnRlcmxvY2tsb3N0XCIsIHBvaW50ZXJMb2NrQ2hhbmdlLCBmYWxzZSlcbiAgICBkb2N1bWVudC5hZGRFdmVudExpc3RlbmVyKFwid2Via2l0cG9pbnRlcmxvY2tsb3N0XCIsIHBvaW50ZXJMb2NrQ2hhbmdlLCBmYWxzZSlcbiAgICBkb2N1bWVudC5hZGRFdmVudExpc3RlbmVyKFwibW96cG9pbnRlcmxvY2tsb3N0XCIsIHBvaW50ZXJMb2NrQ2hhbmdlLCBmYWxzZSlcbiAgICBcbiAgICAvL1VwZGF0ZSBmbGFnc1xuICAgIHNoZWxsLmZ1bGxzY3JlZW4gPSB1c2VGdWxsc2NyZWVuXG4gICAgc2hlbGwucG9pbnRlckxvY2sgPSB1c2VQb2ludGVyTG9ja1xuICBcbiAgICAvL0RlZmF1bHQgbW91c2UgYnV0dG9uIGFsaWFzZXNcbiAgICBzaGVsbC5iaW5kKFwibW91c2UtbGVmdFwiLCAgIFwibW91c2UtMVwiKVxuICAgIHNoZWxsLmJpbmQoXCJtb3VzZS1yaWdodFwiLCAgXCJtb3VzZS0zXCIpXG4gICAgc2hlbGwuYmluZChcIm1vdXNlLW1pZGRsZVwiLCBcIm1vdXNlLTJcIilcbiAgICBcbiAgICAvL0luaXRpYWxpemUgdGljayBjb3VudGVyXG4gICAgc2hlbGwuX2xhc3RUaWNrID0gaHJ0aW1lKClcbiAgICBzaGVsbC5zdGFydFRpbWUgPSBocnRpbWUoKVxuXG4gICAgLy9VbnBhdXNlIHNoZWxsXG4gICAgc2hlbGwucGF1c2VkID0gZmFsc2VcbiAgICBcbiAgICAvL0VtaXQgaW5pdGlhbGl6ZSBldmVudFxuICAgIHNoZWxsLmVtaXQoXCJpbml0XCIpXG4gIH0pfSwgMClcbiAgXG4gIHJldHVybiBzaGVsbFxufVxuXG5tb2R1bGUuZXhwb3J0cyA9IGNyZWF0ZVNoZWxsXG4iLCIvKiAoVGhlIE1JVCBMaWNlbnNlKVxyXG4gKlxyXG4gKiBDb3B5cmlnaHQgKGMpIDIwMTIgQnJhbmRvbiBCZW52aWUgPGh0dHA6Ly9iYmVudmllLmNvbT5cclxuICpcclxuICogUGVybWlzc2lvbiBpcyBoZXJlYnkgZ3JhbnRlZCwgZnJlZSBvZiBjaGFyZ2UsIHRvIGFueSBwZXJzb24gb2J0YWluaW5nIGEgY29weSBvZiB0aGlzIHNvZnR3YXJlIGFuZFxyXG4gKiBhc3NvY2lhdGVkIGRvY3VtZW50YXRpb24gZmlsZXMgKHRoZSAnU29mdHdhcmUnKSwgdG8gZGVhbCBpbiB0aGUgU29mdHdhcmUgd2l0aG91dCByZXN0cmljdGlvbixcclxuICogaW5jbHVkaW5nIHdpdGhvdXQgbGltaXRhdGlvbiB0aGUgcmlnaHRzIHRvIHVzZSwgY29weSwgbW9kaWZ5LCBtZXJnZSwgcHVibGlzaCwgZGlzdHJpYnV0ZSxcclxuICogc3VibGljZW5zZSwgYW5kL29yIHNlbGwgY29waWVzIG9mIHRoZSBTb2Z0d2FyZSwgYW5kIHRvIHBlcm1pdCBwZXJzb25zIHRvIHdob20gdGhlIFNvZnR3YXJlIGlzXHJcbiAqIGZ1cm5pc2hlZCB0byBkbyBzbywgc3ViamVjdCB0byB0aGUgZm9sbG93aW5nIGNvbmRpdGlvbnM6XHJcbiAqXHJcbiAqIFRoZSBhYm92ZSBjb3B5cmlnaHQgbm90aWNlIGFuZCB0aGlzIHBlcm1pc3Npb24gbm90aWNlIHNoYWxsIGJlIGluY2x1ZGVkIHdpdGggYWxsIGNvcGllcyBvclxyXG4gKiBzdWJzdGFudGlhbCBwb3J0aW9ucyBvZiB0aGUgU29mdHdhcmUuXHJcbiAqXHJcbiAqIFRIRSBTT0ZUV0FSRSBJUyBQUk9WSURFRCAnQVMgSVMnLCBXSVRIT1VUIFdBUlJBTlRZIE9GIEFOWSBLSU5ELCBFWFBSRVNTIE9SIElNUExJRUQsIElOQ0xVRElOR1xyXG4gKiBCVVQgTk9UIExJTUlURUQgVE8gVEhFIFdBUlJBTlRJRVMgT0YgTUVSQ0hBTlRBQklMSVRZLCBGSVRORVNTIEZPUiBBIFBBUlRJQ1VMQVIgUFVSUE9TRSBBTkRcclxuICogTk9OSU5GUklOR0VNRU5ULiBJTiBOTyBFVkVOVCBTSEFMTCBUSEUgQVVUSE9SUyBPUiBDT1BZUklHSFQgSE9MREVSUyBCRSBMSUFCTEUgRk9SIEFOWSAgQ0xBSU0sXHJcbiAqIERBTUFHRVMgT1IgT1RIRVIgTElBQklMSVRZLCBXSEVUSEVSIElOIEFOIEFDVElPTiBPRiBDT05UUkFDVCwgVE9SVCBPUiBPVEhFUldJU0UsIEFSSVNJTkcgRlJPTSxcclxuICogT1VUIE9GIE9SIElOIENPTk5FQ1RJT04gV0lUSCBUSEUgU09GVFdBUkUgT1IgVEhFIFVTRSBPUiBPVEhFUiBERUFMSU5HUyBJTiBUSEUgU09GVFdBUkUuXHJcbiAqL1xyXG5cclxuLy8gT3JpZ2luYWwgV2Vha01hcCBpbXBsZW1lbnRhdGlvbiBieSBHb3phbGEgQCBodHRwczovL2dpc3QuZ2l0aHViLmNvbS8xMjY5OTkxXHJcbi8vIFVwZGF0ZWQgYW5kIGJ1Z2ZpeGVkIGJ5IFJheW5vcyBAIGh0dHBzOi8vZ2lzdC5naXRodWIuY29tLzE2MzgwNTlcclxuLy8gRXhwYW5kZWQgYnkgQmVudmllIEAgaHR0cHM6Ly9naXRodWIuY29tL0JlbnZpZS9oYXJtb255LWNvbGxlY3Rpb25zXHJcblxyXG52b2lkIGZ1bmN0aW9uKHN0cmluZ18sIG9iamVjdF8sIGZ1bmN0aW9uXywgcHJvdG90eXBlXywgdG9TdHJpbmdfLFxyXG4gICAgICAgICAgICAgIEFycmF5LCBPYmplY3QsIEZ1bmN0aW9uLCBGUCwgZ2xvYmFsLCBleHBvcnRzLCB1bmRlZmluZWRfLCB1bmRlZmluZWQpe1xyXG5cclxuICB2YXIgZ2V0UHJvcGVydGllcyA9IE9iamVjdC5nZXRPd25Qcm9wZXJ0eU5hbWVzLFxyXG4gICAgICBlczUgPSB0eXBlb2YgZ2V0UHJvcGVydGllcyA9PT0gZnVuY3Rpb25fICYmICEocHJvdG90eXBlXyBpbiBnZXRQcm9wZXJ0aWVzKTtcclxuXHJcbiAgdmFyIGNhbGxiaW5kID0gRlAuYmluZFxyXG4gICAgPyBGUC5iaW5kLmJpbmQoRlAuY2FsbClcclxuICAgIDogKGZ1bmN0aW9uKGNhbGwpe1xyXG4gICAgICAgIHJldHVybiBmdW5jdGlvbihmdW5jKXtcclxuICAgICAgICAgIHJldHVybiBmdW5jdGlvbigpe1xyXG4gICAgICAgICAgICByZXR1cm4gY2FsbC5hcHBseShmdW5jLCBhcmd1bWVudHMpO1xyXG4gICAgICAgICAgfTtcclxuICAgICAgICB9O1xyXG4gICAgICB9KEZQLmNhbGwpKTtcclxuXHJcbiAgdmFyIGZ1bmN0aW9uVG9TdHJpbmcgPSBjYWxsYmluZChGUFt0b1N0cmluZ19dKSxcclxuICAgICAgb2JqZWN0VG9TdHJpbmcgPSBjYWxsYmluZCh7fVt0b1N0cmluZ19dKSxcclxuICAgICAgbnVtYmVyVG9TdHJpbmcgPSBjYWxsYmluZCguMC50b1N0cmluZyksXHJcbiAgICAgIGNhbGwgPSBjYWxsYmluZChGUC5jYWxsKSxcclxuICAgICAgYXBwbHkgPSBjYWxsYmluZChGUC5hcHBseSksXHJcbiAgICAgIGhhc093biA9IGNhbGxiaW5kKHt9Lmhhc093blByb3BlcnR5KSxcclxuICAgICAgcHVzaCA9IGNhbGxiaW5kKFtdLnB1c2gpLFxyXG4gICAgICBzcGxpY2UgPSBjYWxsYmluZChbXS5zcGxpY2UpO1xyXG5cclxuICB2YXIgbmFtZSA9IGZ1bmN0aW9uKGZ1bmMpe1xyXG4gICAgaWYgKHR5cGVvZiBmdW5jICE9PSBmdW5jdGlvbl8pXHJcbiAgICAgIHJldHVybiAnJztcclxuICAgIGVsc2UgaWYgKCduYW1lJyBpbiBmdW5jKVxyXG4gICAgICByZXR1cm4gZnVuYy5uYW1lO1xyXG5cclxuICAgIHJldHVybiBmdW5jdGlvblRvU3RyaW5nKGZ1bmMpLm1hdGNoKC9eXFxuP2Z1bmN0aW9uXFxzPyhcXHcqKT9fP1xcKC8pWzFdO1xyXG4gIH07XHJcblxyXG4gIHZhciBjcmVhdGUgPSBlczVcclxuICAgID8gT2JqZWN0LmNyZWF0ZVxyXG4gICAgOiBmdW5jdGlvbihwcm90bywgZGVzY3Mpe1xyXG4gICAgICAgIHZhciBDdG9yID0gZnVuY3Rpb24oKXt9O1xyXG4gICAgICAgIEN0b3JbcHJvdG90eXBlX10gPSBPYmplY3QocHJvdG8pO1xyXG4gICAgICAgIHZhciBvYmplY3QgPSBuZXcgQ3RvcjtcclxuXHJcbiAgICAgICAgaWYgKGRlc2NzKVxyXG4gICAgICAgICAgZm9yICh2YXIga2V5IGluIGRlc2NzKVxyXG4gICAgICAgICAgICBkZWZpbmVQcm9wZXJ0eShvYmplY3QsIGtleSwgZGVzY3Nba10pO1xyXG5cclxuICAgICAgICByZXR1cm4gb2JqZWN0O1xyXG4gICAgICB9O1xyXG5cclxuXHJcbiAgZnVuY3Rpb24gSGFzaCgpe31cclxuXHJcbiAgaWYgKGVzNSkge1xyXG4gICAgdm9pZCBmdW5jdGlvbihPYmplY3RDcmVhdGUpe1xyXG4gICAgICBIYXNoLnByb3RvdHlwZSA9IE9iamVjdENyZWF0ZShudWxsKTtcclxuICAgICAgZnVuY3Rpb24gaW5oZXJpdChvYmope1xyXG4gICAgICAgIHJldHVybiBPYmplY3RDcmVhdGUob2JqKTtcclxuICAgICAgfVxyXG4gICAgICBIYXNoLmluaGVyaXQgPSBpbmhlcml0O1xyXG4gICAgfShPYmplY3QuY3JlYXRlKTtcclxuICB9IGVsc2Uge1xyXG4gICAgdm9pZCBmdW5jdGlvbihGKXtcclxuICAgICAgdmFyIGlmcmFtZSA9IGRvY3VtZW50LmNyZWF0ZUVsZW1lbnQoJ2lmcmFtZScpO1xyXG4gICAgICBpZnJhbWUuc3R5bGUuZGlzcGxheSA9ICdub25lJztcclxuICAgICAgZG9jdW1lbnQuYm9keS5hcHBlbmRDaGlsZChpZnJhbWUpO1xyXG4gICAgICBpZnJhbWUuc3JjID0gJ2phdmFzY3JpcHQ6J1xyXG4gICAgICBIYXNoLnByb3RvdHlwZSA9IGlmcmFtZS5jb250ZW50V2luZG93Lk9iamVjdC5wcm90b3R5cGU7XHJcbiAgICAgIGRvY3VtZW50LmJvZHkucmVtb3ZlQ2hpbGQoaWZyYW1lKTtcclxuICAgICAgaWZyYW1lID0gbnVsbDtcclxuXHJcbiAgICAgIHZhciBwcm9wcyA9IFsnY29uc3RydWN0b3InLCAnaGFzT3duUHJvcGVydHknLCAncHJvcGVydHlJc0VudW1lcmFibGUnLFxyXG4gICAgICAgICAgICAgICAgICAgJ2lzUHJvdG95cGVPZicsICd0b0xvY2FsZVN0cmluZycsICd0b1N0cmluZycsICd2YWx1ZU9mJ107XHJcblxyXG4gICAgICBmb3IgKHZhciBpPTA7IGkgPCBwcm9wcy5sZW5ndGg7IGkrKylcclxuICAgICAgICBkZWxldGUgSGFzaC5wcm90b3R5cGVbcHJvcHNbaV1dO1xyXG5cclxuICAgICAgZnVuY3Rpb24gaW5oZXJpdChvYmope1xyXG4gICAgICAgIEYucHJvdG90eXBlID0gb2JqO1xyXG4gICAgICAgIG9iaiA9IG5ldyBGO1xyXG4gICAgICAgIEYucHJvdG90eXBlID0gbnVsbDtcclxuICAgICAgICByZXR1cm4gb2JqO1xyXG4gICAgICB9XHJcblxyXG4gICAgICBIYXNoLmluaGVyaXQgPSBpbmhlcml0O1xyXG4gICAgfShmdW5jdGlvbigpe30pO1xyXG4gIH1cclxuXHJcbiAgdmFyIGRlZmluZVByb3BlcnR5ID0gZXM1XHJcbiAgICA/IE9iamVjdC5kZWZpbmVQcm9wZXJ0eVxyXG4gICAgOiBmdW5jdGlvbihvYmplY3QsIGtleSwgZGVzYykge1xyXG4gICAgICAgIG9iamVjdFtrZXldID0gZGVzYy52YWx1ZTtcclxuICAgICAgICByZXR1cm4gb2JqZWN0O1xyXG4gICAgICB9O1xyXG5cclxuICB2YXIgZGVmaW5lID0gZnVuY3Rpb24ob2JqZWN0LCBrZXksIHZhbHVlKXtcclxuICAgIGlmICh0eXBlb2Yga2V5ID09PSBmdW5jdGlvbl8pIHtcclxuICAgICAgdmFsdWUgPSBrZXk7XHJcbiAgICAgIGtleSA9IG5hbWUodmFsdWUpLnJlcGxhY2UoL18kLywgJycpO1xyXG4gICAgfVxyXG5cclxuICAgIHJldHVybiBkZWZpbmVQcm9wZXJ0eShvYmplY3QsIGtleSwgeyBjb25maWd1cmFibGU6IHRydWUsIHdyaXRhYmxlOiB0cnVlLCB2YWx1ZTogdmFsdWUgfSk7XHJcbiAgfTtcclxuXHJcbiAgdmFyIGlzQXJyYXkgPSBlczVcclxuICAgID8gKGZ1bmN0aW9uKGlzQXJyYXkpe1xyXG4gICAgICAgIHJldHVybiBmdW5jdGlvbihvKXtcclxuICAgICAgICAgIHJldHVybiBpc0FycmF5KG8pIHx8IG8gaW5zdGFuY2VvZiBBcnJheTtcclxuICAgICAgICB9O1xyXG4gICAgICB9KShBcnJheS5pc0FycmF5KVxyXG4gICAgOiBmdW5jdGlvbihvKXtcclxuICAgICAgICByZXR1cm4gbyBpbnN0YW5jZW9mIEFycmF5IHx8IG9iamVjdFRvU3RyaW5nKG8pID09PSAnW29iamVjdCBBcnJheV0nO1xyXG4gICAgICB9O1xyXG5cclxuICAvLyAjIyMjIyMjIyMjIyNcclxuICAvLyAjIyMgRGF0YSAjIyNcclxuICAvLyAjIyMjIyMjIyMjIyNcclxuXHJcbiAgdmFyIGJ1aWx0aW5XZWFrTWFwID0gJ1dlYWtNYXAnIGluIGdsb2JhbDtcclxuXHJcbiAgdmFyIE1hcERhdGEgPSBidWlsdGluV2Vha01hcFxyXG4gICAgPyAoZnVuY3Rpb24oKXtcclxuICAgICAgdmFyIEJ1aWx0aW5XZWFrTWFwID0gZ2xvYmFsLldlYWtNYXAsXHJcbiAgICAgICAgICB3bWdldCA9IGNhbGxiaW5kKEJ1aWx0aW5XZWFrTWFwW3Byb3RvdHlwZV9dLmdldCksXHJcbiAgICAgICAgICB3bXNldCA9IGNhbGxiaW5kKEJ1aWx0aW5XZWFrTWFwW3Byb3RvdHlwZV9dLnNldCksXHJcbiAgICAgICAgICB3bWhhcyA9IGNhbGxiaW5kKEJ1aWx0aW5XZWFrTWFwW3Byb3RvdHlwZV9dLmhhcyk7XHJcblxyXG4gICAgICBmdW5jdGlvbiBNYXBEYXRhKG5hbWUpe1xyXG4gICAgICAgIHZhciBtYXAgPSBuZXcgQnVpbHRpbldlYWtNYXA7XHJcblxyXG4gICAgICAgIHRoaXMuZ2V0ID0gZnVuY3Rpb24obyl7XHJcbiAgICAgICAgICByZXR1cm4gd21nZXQobWFwLCBvKTtcclxuICAgICAgICB9O1xyXG4gICAgICAgIHRoaXMuc2V0ID0gZnVuY3Rpb24obywgdil7XHJcbiAgICAgICAgICB3bXNldChtYXAsIG8sIHYpO1xyXG4gICAgICAgIH07XHJcblxyXG4gICAgICAgIGlmIChuYW1lKSB7XHJcbiAgICAgICAgICB0aGlzLndyYXAgPSBmdW5jdGlvbihvLCB2KXtcclxuICAgICAgICAgICAgaWYgKHdtaGFzKG1hcCwgbykpXHJcbiAgICAgICAgICAgICAgdGhyb3cgbmV3IFR5cGVFcnJvcihcIk9iamVjdCBpcyBhbHJlYWR5IGEgXCIgKyBuYW1lKTtcclxuICAgICAgICAgICAgd21zZXQobWFwLCBvLCB2KTtcclxuICAgICAgICAgIH07XHJcbiAgICAgICAgICB0aGlzLnVud3JhcCA9IGZ1bmN0aW9uKG8pe1xyXG4gICAgICAgICAgICB2YXIgc3RvcmFnZSA9IHdtZ2V0KG1hcCwgbyk7XHJcbiAgICAgICAgICAgIGlmICghc3RvcmFnZSlcclxuICAgICAgICAgICAgICB0aHJvdyBuZXcgVHlwZUVycm9yKG5hbWUgKyBcIiBpcyBub3QgZ2VuZXJpY1wiKTtcclxuICAgICAgICAgICAgcmV0dXJuIHN0b3JhZ2U7XHJcbiAgICAgICAgICB9O1xyXG4gICAgICAgIH1cclxuICAgICAgfVxyXG5cclxuICAgICAgcmV0dXJuIE1hcERhdGE7XHJcbiAgICB9KSgpXHJcbiAgICA6IChmdW5jdGlvbigpe1xyXG4gICAgICB2YXIgbG9ja2VyID0gJ3JldHVybiBmdW5jdGlvbihrKXtpZihrPT09cylyZXR1cm4gbH0nLFxyXG4gICAgICAgICAgcmFuZG9tID0gTWF0aC5yYW5kb20sXHJcbiAgICAgICAgICB1aWRzID0gbmV3IEhhc2gsXHJcbiAgICAgICAgICBzbGljZSA9IGNhbGxiaW5kKCcnLnNsaWNlKSxcclxuICAgICAgICAgIGluZGV4T2YgPSBjYWxsYmluZChbXS5pbmRleE9mKTtcclxuXHJcbiAgICAgIHZhciBjcmVhdGVVSUQgPSBmdW5jdGlvbigpe1xyXG4gICAgICAgIHZhciBrZXkgPSBzbGljZShudW1iZXJUb1N0cmluZyhyYW5kb20oKSwgMzYpLCAyKTtcclxuICAgICAgICByZXR1cm4ga2V5IGluIHVpZHMgPyBjcmVhdGVVSUQoKSA6IHVpZHNba2V5XSA9IGtleTtcclxuICAgICAgfTtcclxuXHJcbiAgICAgIHZhciBnbG9iYWxJRCA9IGNyZWF0ZVVJRCgpO1xyXG5cclxuICAgICAgLy8gY29tbW9uIHBlci1vYmplY3Qgc3RvcmFnZSBhcmVhIG1hZGUgdmlzaWJsZSBieSBwYXRjaGluZyBnZXRPd25Qcm9wZXJ0eU5hbWVzJ1xyXG4gICAgICBmdW5jdGlvbiBnZXRPd25Qcm9wZXJ0eU5hbWVzKG9iail7XHJcbiAgICAgICAgdmFyIHByb3BzID0gZ2V0UHJvcGVydGllcyhvYmopO1xyXG4gICAgICAgIGlmIChoYXNPd24ob2JqLCBnbG9iYWxJRCkpXHJcbiAgICAgICAgICBzcGxpY2UocHJvcHMsIGluZGV4T2YocHJvcHMsIGdsb2JhbElEKSwgMSk7XHJcbiAgICAgICAgcmV0dXJuIHByb3BzO1xyXG4gICAgICB9XHJcblxyXG4gICAgICBpZiAoZXM1KSB7XHJcbiAgICAgICAgLy8gY2hlY2sgZm9yIHRoZSByYW5kb20ga2V5IG9uIGFuIG9iamVjdCwgY3JlYXRlIG5ldyBzdG9yYWdlIGlmIG1pc3NpbmcsIHJldHVybiBpdFxyXG4gICAgICAgIHZhciBzdG9yYWdlID0gZnVuY3Rpb24ob2JqKXtcclxuICAgICAgICAgIGlmICghaGFzT3duKG9iaiwgZ2xvYmFsSUQpKVxyXG4gICAgICAgICAgICBkZWZpbmVQcm9wZXJ0eShvYmosIGdsb2JhbElELCB7IHZhbHVlOiBuZXcgSGFzaCB9KTtcclxuICAgICAgICAgIHJldHVybiBvYmpbZ2xvYmFsSURdO1xyXG4gICAgICAgIH07XHJcblxyXG4gICAgICAgIGRlZmluZShPYmplY3QsIGdldE93blByb3BlcnR5TmFtZXMpO1xyXG4gICAgICB9IGVsc2Uge1xyXG5cclxuICAgICAgICB2YXIgdG9TdHJpbmdUb1N0cmluZyA9IGZ1bmN0aW9uKHMpe1xyXG4gICAgICAgICAgZnVuY3Rpb24gdG9TdHJpbmcoKXsgcmV0dXJuIHMgfVxyXG4gICAgICAgICAgcmV0dXJuIHRvU3RyaW5nW3RvU3RyaW5nX10gPSB0b1N0cmluZztcclxuICAgICAgICB9KE9iamVjdFtwcm90b3R5cGVfXVt0b1N0cmluZ19dKycnKTtcclxuXHJcbiAgICAgICAgLy8gc3RvcmUgdGhlIHZhbHVlcyBvbiBhIGN1c3RvbSB2YWx1ZU9mIGluIG9yZGVyIHRvIGhpZGUgdGhlbSBidXQgc3RvcmUgdGhlbSBsb2NhbGx5XHJcbiAgICAgICAgdmFyIHN0b3JhZ2UgPSBmdW5jdGlvbihvYmope1xyXG4gICAgICAgICAgaWYgKGhhc093bihvYmosIHRvU3RyaW5nXykgJiYgZ2xvYmFsSUQgaW4gb2JqW3RvU3RyaW5nX10pXHJcbiAgICAgICAgICAgIHJldHVybiBvYmpbdG9TdHJpbmdfXVtnbG9iYWxJRF07XHJcblxyXG4gICAgICAgICAgaWYgKCEodG9TdHJpbmdfIGluIG9iaikpXHJcbiAgICAgICAgICAgIHRocm93IG5ldyBFcnJvcihcIkNhbid0IHN0b3JlIHZhbHVlcyBmb3IgXCIrb2JqKTtcclxuXHJcbiAgICAgICAgICB2YXIgb2xkVG9TdHJpbmcgPSBvYmpbdG9TdHJpbmdfXTtcclxuICAgICAgICAgIGZ1bmN0aW9uIHRvU3RyaW5nKCl7IHJldHVybiBvbGRUb1N0cmluZy5jYWxsKHRoaXMpIH1cclxuICAgICAgICAgIG9ialt0b1N0cmluZ19dID0gdG9TdHJpbmc7XHJcbiAgICAgICAgICB0b1N0cmluZ1t0b1N0cmluZ19dID0gdG9TdHJpbmdUb1N0cmluZztcclxuICAgICAgICAgIHJldHVybiB0b1N0cmluZ1tnbG9iYWxJRF0gPSB7fTtcclxuICAgICAgICB9O1xyXG4gICAgICB9XHJcblxyXG5cclxuXHJcbiAgICAgIC8vIHNoaW0gZm9yIFtbTWFwRGF0YV1dIGZyb20gZXM2IHNwZWMsIGFuZCBwdWxscyBkb3VibGUgZHV0eSBhcyBXZWFrTWFwIHN0b3JhZ2VcclxuICAgICAgZnVuY3Rpb24gTWFwRGF0YShuYW1lKXtcclxuICAgICAgICB2YXIgcHVpZCA9IGNyZWF0ZVVJRCgpLFxyXG4gICAgICAgICAgICBpdWlkID0gY3JlYXRlVUlEKCksXHJcbiAgICAgICAgICAgIHNlY3JldCA9IHsgdmFsdWU6IHVuZGVmaW5lZCB9O1xyXG5cclxuICAgICAgICB2YXIgYXR0YWNoID0gZnVuY3Rpb24ob2JqKXtcclxuICAgICAgICAgIHZhciBzdG9yZSA9IHN0b3JhZ2Uob2JqKTtcclxuICAgICAgICAgIGlmIChoYXNPd24oc3RvcmUsIHB1aWQpKVxyXG4gICAgICAgICAgICByZXR1cm4gc3RvcmVbcHVpZF0oc2VjcmV0KTtcclxuXHJcbiAgICAgICAgICB2YXIgbG9ja2JveCA9IG5ldyBIYXNoO1xyXG4gICAgICAgICAgZGVmaW5lUHJvcGVydHkobG9ja2JveCwgaXVpZCwgc2VjcmV0KTtcclxuICAgICAgICAgIGRlZmluZVByb3BlcnR5KHN0b3JlLCBwdWlkLCB7XHJcbiAgICAgICAgICAgIHZhbHVlOiBuZXcgRnVuY3Rpb24oJ3MnLCAnbCcsIGxvY2tlcikoc2VjcmV0LCBsb2NrYm94KVxyXG4gICAgICAgICAgfSk7XHJcbiAgICAgICAgICByZXR1cm4gbG9ja2JveDtcclxuICAgICAgICB9O1xyXG5cclxuICAgICAgICB0aGlzLmdldCA9IGZ1bmN0aW9uKG8pe1xyXG4gICAgICAgICAgcmV0dXJuIGF0dGFjaChvKVtpdWlkXTtcclxuICAgICAgICB9O1xyXG4gICAgICAgIHRoaXMuc2V0ID0gZnVuY3Rpb24obywgdil7XHJcbiAgICAgICAgICBhdHRhY2gobylbaXVpZF0gPSB2O1xyXG4gICAgICAgIH07XHJcblxyXG4gICAgICAgIGlmIChuYW1lKSB7XHJcbiAgICAgICAgICB0aGlzLndyYXAgPSBmdW5jdGlvbihvLCB2KXtcclxuICAgICAgICAgICAgdmFyIGxvY2tib3ggPSBhdHRhY2gobyk7XHJcbiAgICAgICAgICAgIGlmIChsb2NrYm94W2l1aWRdKVxyXG4gICAgICAgICAgICAgIHRocm93IG5ldyBUeXBlRXJyb3IoXCJPYmplY3QgaXMgYWxyZWFkeSBhIFwiICsgbmFtZSk7XHJcbiAgICAgICAgICAgIGxvY2tib3hbaXVpZF0gPSB2O1xyXG4gICAgICAgICAgfTtcclxuICAgICAgICAgIHRoaXMudW53cmFwID0gZnVuY3Rpb24obyl7XHJcbiAgICAgICAgICAgIHZhciBzdG9yYWdlID0gYXR0YWNoKG8pW2l1aWRdO1xyXG4gICAgICAgICAgICBpZiAoIXN0b3JhZ2UpXHJcbiAgICAgICAgICAgICAgdGhyb3cgbmV3IFR5cGVFcnJvcihuYW1lICsgXCIgaXMgbm90IGdlbmVyaWNcIik7XHJcbiAgICAgICAgICAgIHJldHVybiBzdG9yYWdlO1xyXG4gICAgICAgICAgfTtcclxuICAgICAgICB9XHJcbiAgICAgIH1cclxuXHJcbiAgICAgIHJldHVybiBNYXBEYXRhO1xyXG4gICAgfSgpKTtcclxuXHJcbiAgdmFyIGV4cG9ydGVyID0gKGZ1bmN0aW9uKCl7XHJcbiAgICAvLyBbbmF0aXZlIGNvZGVdIGxvb2tzIHNsaWdodGx5IGRpZmZlcmVudCBpbiBlYWNoIGVuZ2luZVxyXG4gICAgdmFyIHNyYyA9ICgnJytPYmplY3QpLnNwbGl0KCdPYmplY3QnKTtcclxuXHJcbiAgICAvLyBmYWtlIFtuYXRpdmUgY29kZV1cclxuICAgIGZ1bmN0aW9uIHRvU3RyaW5nKCl7XHJcbiAgICAgIHJldHVybiBzcmNbMF0gKyBuYW1lKHRoaXMpICsgc3JjWzFdO1xyXG4gICAgfVxyXG5cclxuICAgIGRlZmluZSh0b1N0cmluZywgdG9TdHJpbmcpO1xyXG5cclxuICAgIC8vIGF0dGVtcHQgdG8gdXNlIF9fcHJvdG9fXyBzbyB0aGUgbWV0aG9kcyBkb24ndCBhbGwgaGF2ZSBhbiBvd24gdG9TdHJpbmdcclxuICAgIHZhciBwcmVwRnVuY3Rpb24gPSB7IF9fcHJvdG9fXzogW10gfSBpbnN0YW5jZW9mIEFycmF5XHJcbiAgICAgID8gZnVuY3Rpb24oZnVuYyl7IGZ1bmMuX19wcm90b19fID0gdG9TdHJpbmcgfVxyXG4gICAgICA6IGZ1bmN0aW9uKGZ1bmMpeyBkZWZpbmUoZnVuYywgdG9TdHJpbmcpIH07XHJcblxyXG4gICAgLy8gYXNzZW1ibGUgYW4gYXJyYXkgb2YgZnVuY3Rpb25zIGludG8gYSBmdWxseSBmb3JtZWQgY2xhc3NcclxuICAgIHZhciBwcmVwYXJlID0gZnVuY3Rpb24obWV0aG9kcyl7XHJcbiAgICAgIHZhciBDdG9yID0gbWV0aG9kcy5zaGlmdCgpLFxyXG4gICAgICAgICAgYnJhbmQgPSAnW29iamVjdCAnICsgbmFtZShDdG9yKSArICddJztcclxuXHJcbiAgICAgIGZ1bmN0aW9uIHRvU3RyaW5nKCl7IHJldHVybiBicmFuZCB9XHJcbiAgICAgIG1ldGhvZHMucHVzaCh0b1N0cmluZyk7XHJcbiAgICAgIHByZXBGdW5jdGlvbihDdG9yKTtcclxuXHJcbiAgICAgIGZvciAodmFyIGk9MDsgaSA8IG1ldGhvZHMubGVuZ3RoOyBpKyspIHtcclxuICAgICAgICBwcmVwRnVuY3Rpb24obWV0aG9kc1tpXSk7XHJcbiAgICAgICAgZGVmaW5lKEN0b3JbcHJvdG90eXBlX10sIG1ldGhvZHNbaV0pO1xyXG4gICAgICB9XHJcblxyXG4gICAgICByZXR1cm4gQ3RvcjtcclxuICAgIH07XHJcblxyXG4gICAgcmV0dXJuIGZ1bmN0aW9uKG5hbWUsIGluaXQpe1xyXG4gICAgICBpZiAobmFtZSBpbiBleHBvcnRzKVxyXG4gICAgICAgIHJldHVybiBleHBvcnRzW25hbWVdO1xyXG5cclxuICAgICAgdmFyIGRhdGEgPSBuZXcgTWFwRGF0YShuYW1lKTtcclxuXHJcbiAgICAgIHJldHVybiBleHBvcnRzW25hbWVdID0gcHJlcGFyZShpbml0KFxyXG4gICAgICAgIGZ1bmN0aW9uKGNvbGxlY3Rpb24sIHZhbHVlKXtcclxuICAgICAgICAgIGRhdGEud3JhcChjb2xsZWN0aW9uLCB2YWx1ZSk7XHJcbiAgICAgICAgfSxcclxuICAgICAgICBmdW5jdGlvbihjb2xsZWN0aW9uKXtcclxuICAgICAgICAgIHJldHVybiBkYXRhLnVud3JhcChjb2xsZWN0aW9uKTtcclxuICAgICAgICB9XHJcbiAgICAgICkpO1xyXG4gICAgfTtcclxuICB9KCkpO1xyXG5cclxuXHJcbiAgLy8gaW5pdGlhbGl6ZSBjb2xsZWN0aW9uIHdpdGggYW4gaXRlcmFibGUsIGN1cnJlbnRseSBvbmx5IHN1cHBvcnRzIGZvckVhY2ggZnVuY3Rpb25cclxuICB2YXIgaW5pdGlhbGl6ZSA9IGZ1bmN0aW9uKGl0ZXJhYmxlLCBjYWxsYmFjayl7XHJcbiAgICBpZiAoaXRlcmFibGUgIT09IG51bGwgJiYgdHlwZW9mIGl0ZXJhYmxlID09PSBvYmplY3RfICYmIHR5cGVvZiBpdGVyYWJsZS5mb3JFYWNoID09PSBmdW5jdGlvbl8pIHtcclxuICAgICAgaXRlcmFibGUuZm9yRWFjaChmdW5jdGlvbihpdGVtLCBpKXtcclxuICAgICAgICBpZiAoaXNBcnJheShpdGVtKSAmJiBpdGVtLmxlbmd0aCA9PT0gMilcclxuICAgICAgICAgIGNhbGxiYWNrKGl0ZXJhYmxlW2ldWzBdLCBpdGVyYWJsZVtpXVsxXSk7XHJcbiAgICAgICAgZWxzZVxyXG4gICAgICAgICAgY2FsbGJhY2soaXRlcmFibGVbaV0sIGkpO1xyXG4gICAgICB9KTtcclxuICAgIH1cclxuICB9XHJcblxyXG4gIC8vIGF0dGVtcHQgdG8gZml4IHRoZSBuYW1lIG9mIFwiZGVsZXRlX1wiIG1ldGhvZHMsIHNob3VsZCB3b3JrIGluIFY4IGFuZCBzcGlkZXJtb25rZXlcclxuICB2YXIgZml4RGVsZXRlID0gZnVuY3Rpb24oZnVuYywgc2NvcGVOYW1lcywgc2NvcGVWYWx1ZXMpe1xyXG4gICAgdHJ5IHtcclxuICAgICAgc2NvcGVOYW1lc1tzY29wZU5hbWVzLmxlbmd0aF0gPSAoJ3JldHVybiAnK2Z1bmMpLnJlcGxhY2UoJ2VfJywgJ1xcXFx1MDA2NScpO1xyXG4gICAgICByZXR1cm4gRnVuY3Rpb24uYXBwbHkoMCwgc2NvcGVOYW1lcykuYXBwbHkoMCwgc2NvcGVWYWx1ZXMpO1xyXG4gICAgfSBjYXRjaCAoZSkge1xyXG4gICAgICByZXR1cm4gZnVuYztcclxuICAgIH1cclxuICB9XHJcblxyXG4gIHZhciBXTSwgSE0sIE07XHJcblxyXG4gIC8vICMjIyMjIyMjIyMjIyMjI1xyXG4gIC8vICMjIyBXZWFrTWFwICMjI1xyXG4gIC8vICMjIyMjIyMjIyMjIyMjI1xyXG5cclxuICBXTSA9IGJ1aWx0aW5XZWFrTWFwID8gKGV4cG9ydHMuV2Vha01hcCA9IGdsb2JhbC5XZWFrTWFwKSA6IGV4cG9ydGVyKCdXZWFrTWFwJywgZnVuY3Rpb24od3JhcCwgdW53cmFwKXtcclxuICAgIHZhciBwcm90b3R5cGUgPSBXZWFrTWFwW3Byb3RvdHlwZV9dO1xyXG4gICAgdmFyIHZhbGlkYXRlID0gZnVuY3Rpb24oa2V5KXtcclxuICAgICAgaWYgKGtleSA9PSBudWxsIHx8IHR5cGVvZiBrZXkgIT09IG9iamVjdF8gJiYgdHlwZW9mIGtleSAhPT0gZnVuY3Rpb25fKVxyXG4gICAgICAgIHRocm93IG5ldyBUeXBlRXJyb3IoXCJJbnZhbGlkIFdlYWtNYXAga2V5XCIpO1xyXG4gICAgfTtcclxuXHJcbiAgICAvKipcclxuICAgICAqIEBjbGFzcyAgICAgICAgV2Vha01hcFxyXG4gICAgICogQGRlc2NyaXB0aW9uICBDb2xsZWN0aW9uIHVzaW5nIG9iamVjdHMgd2l0aCB1bmlxdWUgaWRlbnRpdGllcyBhcyBrZXlzIHRoYXQgZGlzYWxsb3dzIGVudW1lcmF0aW9uXHJcbiAgICAgKiAgICAgICAgICAgICAgIGFuZCBhbGxvd3MgZm9yIGJldHRlciBnYXJiYWdlIGNvbGxlY3Rpb24uXHJcbiAgICAgKiBAcGFyYW0gICAgICAgIHtJdGVyYWJsZX0gW2l0ZXJhYmxlXSAgQW4gaXRlbSB0byBwb3B1bGF0ZSB0aGUgY29sbGVjdGlvbiB3aXRoLlxyXG4gICAgICovXHJcbiAgICBmdW5jdGlvbiBXZWFrTWFwKGl0ZXJhYmxlKXtcclxuICAgICAgaWYgKHRoaXMgPT09IGdsb2JhbCB8fCB0aGlzID09IG51bGwgfHwgdGhpcyA9PT0gcHJvdG90eXBlKVxyXG4gICAgICAgIHJldHVybiBuZXcgV2Vha01hcChpdGVyYWJsZSk7XHJcblxyXG4gICAgICB3cmFwKHRoaXMsIG5ldyBNYXBEYXRhKTtcclxuXHJcbiAgICAgIHZhciBzZWxmID0gdGhpcztcclxuICAgICAgaXRlcmFibGUgJiYgaW5pdGlhbGl6ZShpdGVyYWJsZSwgZnVuY3Rpb24odmFsdWUsIGtleSl7XHJcbiAgICAgICAgY2FsbChzZXQsIHNlbGYsIHZhbHVlLCBrZXkpO1xyXG4gICAgICB9KTtcclxuICAgIH1cclxuICAgIC8qKlxyXG4gICAgICogQG1ldGhvZCAgICAgICA8Z2V0PlxyXG4gICAgICogQGRlc2NyaXB0aW9uICBSZXRyaWV2ZSB0aGUgdmFsdWUgaW4gdGhlIGNvbGxlY3Rpb24gdGhhdCBtYXRjaGVzIGtleVxyXG4gICAgICogQHBhcmFtICAgICAgICB7QW55fSBrZXlcclxuICAgICAqIEByZXR1cm4gICAgICAge0FueX1cclxuICAgICAqL1xyXG4gICAgZnVuY3Rpb24gZ2V0KGtleSl7XHJcbiAgICAgIHZhbGlkYXRlKGtleSk7XHJcbiAgICAgIHZhciB2YWx1ZSA9IHVud3JhcCh0aGlzKS5nZXQoa2V5KTtcclxuICAgICAgcmV0dXJuIHZhbHVlID09PSB1bmRlZmluZWRfID8gdW5kZWZpbmVkIDogdmFsdWU7XHJcbiAgICB9XHJcbiAgICAvKipcclxuICAgICAqIEBtZXRob2QgICAgICAgPHNldD5cclxuICAgICAqIEBkZXNjcmlwdGlvbiAgQWRkIG9yIHVwZGF0ZSBhIHBhaXIgaW4gdGhlIGNvbGxlY3Rpb24uIEVuZm9yY2VzIHVuaXF1ZW5lc3MgYnkgb3ZlcndyaXRpbmcuXHJcbiAgICAgKiBAcGFyYW0gICAgICAgIHtBbnl9IGtleVxyXG4gICAgICogQHBhcmFtICAgICAgICB7QW55fSB2YWxcclxuICAgICAqKi9cclxuICAgIGZ1bmN0aW9uIHNldChrZXksIHZhbHVlKXtcclxuICAgICAgdmFsaWRhdGUoa2V5KTtcclxuICAgICAgLy8gc3RvcmUgYSB0b2tlbiBmb3IgZXhwbGljaXQgdW5kZWZpbmVkIHNvIHRoYXQgXCJoYXNcIiB3b3JrcyBjb3JyZWN0bHlcclxuICAgICAgdW53cmFwKHRoaXMpLnNldChrZXksIHZhbHVlID09PSB1bmRlZmluZWQgPyB1bmRlZmluZWRfIDogdmFsdWUpO1xyXG4gICAgfVxyXG4gICAgLypcclxuICAgICAqIEBtZXRob2QgICAgICAgPGhhcz5cclxuICAgICAqIEBkZXNjcmlwdGlvbiAgQ2hlY2sgaWYga2V5IGlzIGluIHRoZSBjb2xsZWN0aW9uXHJcbiAgICAgKiBAcGFyYW0gICAgICAgIHtBbnl9IGtleVxyXG4gICAgICogQHJldHVybiAgICAgICB7Qm9vbGVhbn1cclxuICAgICAqKi9cclxuICAgIGZ1bmN0aW9uIGhhcyhrZXkpe1xyXG4gICAgICB2YWxpZGF0ZShrZXkpO1xyXG4gICAgICByZXR1cm4gdW53cmFwKHRoaXMpLmdldChrZXkpICE9PSB1bmRlZmluZWQ7XHJcbiAgICB9XHJcbiAgICAvKipcclxuICAgICAqIEBtZXRob2QgICAgICAgPGRlbGV0ZT5cclxuICAgICAqIEBkZXNjcmlwdGlvbiAgUmVtb3ZlIGtleSBhbmQgbWF0Y2hpbmcgdmFsdWUgaWYgZm91bmRcclxuICAgICAqIEBwYXJhbSAgICAgICAge0FueX0ga2V5XHJcbiAgICAgKiBAcmV0dXJuICAgICAgIHtCb29sZWFufSB0cnVlIGlmIGl0ZW0gd2FzIGluIGNvbGxlY3Rpb25cclxuICAgICAqL1xyXG4gICAgZnVuY3Rpb24gZGVsZXRlXyhrZXkpe1xyXG4gICAgICB2YWxpZGF0ZShrZXkpO1xyXG4gICAgICB2YXIgZGF0YSA9IHVud3JhcCh0aGlzKTtcclxuXHJcbiAgICAgIGlmIChkYXRhLmdldChrZXkpID09PSB1bmRlZmluZWQpXHJcbiAgICAgICAgcmV0dXJuIGZhbHNlO1xyXG5cclxuICAgICAgZGF0YS5zZXQoa2V5LCB1bmRlZmluZWQpO1xyXG4gICAgICByZXR1cm4gdHJ1ZTtcclxuICAgIH1cclxuXHJcbiAgICBkZWxldGVfID0gZml4RGVsZXRlKGRlbGV0ZV8sIFsndmFsaWRhdGUnLCAndW53cmFwJ10sIFt2YWxpZGF0ZSwgdW53cmFwXSk7XHJcbiAgICByZXR1cm4gW1dlYWtNYXAsIGdldCwgc2V0LCBoYXMsIGRlbGV0ZV9dO1xyXG4gIH0pO1xyXG5cclxuXHJcbiAgLy8gIyMjIyMjIyMjIyMjIyMjXHJcbiAgLy8gIyMjIEhhc2hNYXAgIyMjXHJcbiAgLy8gIyMjIyMjIyMjIyMjIyMjXHJcblxyXG4gIEhNID0gZXhwb3J0ZXIoJ0hhc2hNYXAnLCBmdW5jdGlvbih3cmFwLCB1bndyYXApe1xyXG4gICAgLy8gc2VwYXJhdGUgbnVtYmVycywgc3RyaW5ncywgYW5kIGF0b21zIHRvIGNvbXBlbnNhdGUgZm9yIGtleSBjb2VyY2lvbiB0byBzdHJpbmdcclxuXHJcbiAgICB2YXIgcHJvdG90eXBlID0gSGFzaE1hcFtwcm90b3R5cGVfXSxcclxuICAgICAgICBTVFJJTkcgPSAwLCBOVU1CRVIgPSAxLCBPVEhFUiA9IDIsXHJcbiAgICAgICAgb3RoZXJzID0geyAndHJ1ZSc6IHRydWUsICdmYWxzZSc6IGZhbHNlLCAnbnVsbCc6IG51bGwsIDA6IC0wIH07XHJcblxyXG4gICAgdmFyIHByb3RvID0gTWF0aC5yYW5kb20oKS50b1N0cmluZygzNikuc2xpY2UoMik7XHJcblxyXG4gICAgdmFyIGNvZXJjZSA9IGZ1bmN0aW9uKGtleSl7XHJcbiAgICAgIHJldHVybiBrZXkgPT09ICdfX3Byb3RvX18nID8gcHJvdG8gOiBrZXk7XHJcbiAgICB9O1xyXG5cclxuICAgIHZhciB1bmNvZXJjZSA9IGZ1bmN0aW9uKHR5cGUsIGtleSl7XHJcbiAgICAgIHN3aXRjaCAodHlwZSkge1xyXG4gICAgICAgIGNhc2UgU1RSSU5HOiByZXR1cm4ga2V5ID09PSBwcm90byA/ICdfX3Byb3RvX18nIDoga2V5O1xyXG4gICAgICAgIGNhc2UgTlVNQkVSOiByZXR1cm4gK2tleTtcclxuICAgICAgICBjYXNlIE9USEVSOiByZXR1cm4gb3RoZXJzW2tleV07XHJcbiAgICAgIH1cclxuICAgIH1cclxuXHJcblxyXG4gICAgdmFyIHZhbGlkYXRlID0gZnVuY3Rpb24oa2V5KXtcclxuICAgICAgaWYgKGtleSA9PSBudWxsKSByZXR1cm4gT1RIRVI7XHJcbiAgICAgIHN3aXRjaCAodHlwZW9mIGtleSkge1xyXG4gICAgICAgIGNhc2UgJ2Jvb2xlYW4nOiByZXR1cm4gT1RIRVI7XHJcbiAgICAgICAgY2FzZSBzdHJpbmdfOiByZXR1cm4gU1RSSU5HO1xyXG4gICAgICAgIC8vIG5lZ2F0aXZlIHplcm8gaGFzIHRvIGJlIGV4cGxpY2l0bHkgYWNjb3VudGVkIGZvclxyXG4gICAgICAgIGNhc2UgJ251bWJlcic6IHJldHVybiBrZXkgPT09IDAgJiYgSW5maW5pdHkgLyBrZXkgPT09IC1JbmZpbml0eSA/IE9USEVSIDogTlVNQkVSO1xyXG4gICAgICAgIGRlZmF1bHQ6IHRocm93IG5ldyBUeXBlRXJyb3IoXCJJbnZhbGlkIEhhc2hNYXAga2V5XCIpO1xyXG4gICAgICB9XHJcbiAgICB9XHJcblxyXG4gICAgLyoqXHJcbiAgICAgKiBAY2xhc3MgICAgICAgICAgSGFzaE1hcFxyXG4gICAgICogQGRlc2NyaXB0aW9uICAgIENvbGxlY3Rpb24gdGhhdCBvbmx5IGFsbG93cyBwcmltaXRpdmVzIHRvIGJlIGtleXMuXHJcbiAgICAgKiBAcGFyYW0gICAgICAgICAge0l0ZXJhYmxlfSBbaXRlcmFibGVdICBBbiBpdGVtIHRvIHBvcHVsYXRlIHRoZSBjb2xsZWN0aW9uIHdpdGguXHJcbiAgICAgKi9cclxuICAgIGZ1bmN0aW9uIEhhc2hNYXAoaXRlcmFibGUpe1xyXG4gICAgICBpZiAodGhpcyA9PT0gZ2xvYmFsIHx8IHRoaXMgPT0gbnVsbCB8fCB0aGlzID09PSBwcm90b3R5cGUpXHJcbiAgICAgICAgcmV0dXJuIG5ldyBIYXNoTWFwKGl0ZXJhYmxlKTtcclxuXHJcbiAgICAgIHdyYXAodGhpcywge1xyXG4gICAgICAgIHNpemU6IDAsXHJcbiAgICAgICAgMDogbmV3IEhhc2gsXHJcbiAgICAgICAgMTogbmV3IEhhc2gsXHJcbiAgICAgICAgMjogbmV3IEhhc2hcclxuICAgICAgfSk7XHJcblxyXG4gICAgICB2YXIgc2VsZiA9IHRoaXM7XHJcbiAgICAgIGl0ZXJhYmxlICYmIGluaXRpYWxpemUoaXRlcmFibGUsIGZ1bmN0aW9uKHZhbHVlLCBrZXkpe1xyXG4gICAgICAgIGNhbGwoc2V0LCBzZWxmLCB2YWx1ZSwga2V5KTtcclxuICAgICAgfSk7XHJcbiAgICB9XHJcbiAgICAvKipcclxuICAgICAqIEBtZXRob2QgICAgICAgPGdldD5cclxuICAgICAqIEBkZXNjcmlwdGlvbiAgUmV0cmlldmUgdGhlIHZhbHVlIGluIHRoZSBjb2xsZWN0aW9uIHRoYXQgbWF0Y2hlcyBrZXlcclxuICAgICAqIEBwYXJhbSAgICAgICAge0FueX0ga2V5XHJcbiAgICAgKiBAcmV0dXJuICAgICAgIHtBbnl9XHJcbiAgICAgKi9cclxuICAgIGZ1bmN0aW9uIGdldChrZXkpe1xyXG4gICAgICByZXR1cm4gdW53cmFwKHRoaXMpW3ZhbGlkYXRlKGtleSldW2NvZXJjZShrZXkpXTtcclxuICAgIH1cclxuICAgIC8qKlxyXG4gICAgICogQG1ldGhvZCAgICAgICA8c2V0PlxyXG4gICAgICogQGRlc2NyaXB0aW9uICBBZGQgb3IgdXBkYXRlIGEgcGFpciBpbiB0aGUgY29sbGVjdGlvbi4gRW5mb3JjZXMgdW5pcXVlbmVzcyBieSBvdmVyd3JpdGluZy5cclxuICAgICAqIEBwYXJhbSAgICAgICAge0FueX0ga2V5XHJcbiAgICAgKiBAcGFyYW0gICAgICAgIHtBbnl9IHZhbFxyXG4gICAgICoqL1xyXG4gICAgZnVuY3Rpb24gc2V0KGtleSwgdmFsdWUpe1xyXG4gICAgICB2YXIgaXRlbXMgPSB1bndyYXAodGhpcyksXHJcbiAgICAgICAgICBkYXRhID0gaXRlbXNbdmFsaWRhdGUoa2V5KV07XHJcblxyXG4gICAgICBrZXkgPSBjb2VyY2Uoa2V5KTtcclxuICAgICAga2V5IGluIGRhdGEgfHwgaXRlbXMuc2l6ZSsrO1xyXG4gICAgICBkYXRhW2tleV0gPSB2YWx1ZTtcclxuICAgIH1cclxuICAgIC8qKlxyXG4gICAgICogQG1ldGhvZCAgICAgICA8aGFzPlxyXG4gICAgICogQGRlc2NyaXB0aW9uICBDaGVjayBpZiBrZXkgZXhpc3RzIGluIHRoZSBjb2xsZWN0aW9uLlxyXG4gICAgICogQHBhcmFtICAgICAgICB7QW55fSBrZXlcclxuICAgICAqIEByZXR1cm4gICAgICAge0Jvb2xlYW59IGlzIGluIGNvbGxlY3Rpb25cclxuICAgICAqKi9cclxuICAgIGZ1bmN0aW9uIGhhcyhrZXkpe1xyXG4gICAgICByZXR1cm4gY29lcmNlKGtleSkgaW4gdW53cmFwKHRoaXMpW3ZhbGlkYXRlKGtleSldO1xyXG4gICAgfVxyXG4gICAgLyoqXHJcbiAgICAgKiBAbWV0aG9kICAgICAgIDxkZWxldGU+XHJcbiAgICAgKiBAZGVzY3JpcHRpb24gIFJlbW92ZSBrZXkgYW5kIG1hdGNoaW5nIHZhbHVlIGlmIGZvdW5kXHJcbiAgICAgKiBAcGFyYW0gICAgICAgIHtBbnl9IGtleVxyXG4gICAgICogQHJldHVybiAgICAgICB7Qm9vbGVhbn0gdHJ1ZSBpZiBpdGVtIHdhcyBpbiBjb2xsZWN0aW9uXHJcbiAgICAgKi9cclxuICAgIGZ1bmN0aW9uIGRlbGV0ZV8oa2V5KXtcclxuICAgICAgdmFyIGl0ZW1zID0gdW53cmFwKHRoaXMpLFxyXG4gICAgICAgICAgZGF0YSA9IGl0ZW1zW3ZhbGlkYXRlKGtleSldO1xyXG5cclxuICAgICAga2V5ID0gY29lcmNlKGtleSk7XHJcbiAgICAgIGlmIChrZXkgaW4gZGF0YSkge1xyXG4gICAgICAgIGRlbGV0ZSBkYXRhW2tleV07XHJcbiAgICAgICAgaXRlbXMuc2l6ZS0tO1xyXG4gICAgICAgIHJldHVybiB0cnVlO1xyXG4gICAgICB9XHJcblxyXG4gICAgICByZXR1cm4gZmFsc2U7XHJcbiAgICB9XHJcbiAgICAvKipcclxuICAgICAqIEBtZXRob2QgICAgICAgPHNpemU+XHJcbiAgICAgKiBAZGVzY3JpcHRpb24gIFJldHJpZXZlIHRoZSBhbW91bnQgb2YgaXRlbXMgaW4gdGhlIGNvbGxlY3Rpb25cclxuICAgICAqIEByZXR1cm4gICAgICAge051bWJlcn1cclxuICAgICAqL1xyXG4gICAgZnVuY3Rpb24gc2l6ZSgpe1xyXG4gICAgICByZXR1cm4gdW53cmFwKHRoaXMpLnNpemU7XHJcbiAgICB9XHJcbiAgICAvKipcclxuICAgICAqIEBtZXRob2QgICAgICAgPGZvckVhY2g+XHJcbiAgICAgKiBAZGVzY3JpcHRpb24gIExvb3AgdGhyb3VnaCB0aGUgY29sbGVjdGlvbiByYWlzaW5nIGNhbGxiYWNrIGZvciBlYWNoXHJcbiAgICAgKiBAcGFyYW0gICAgICAgIHtGdW5jdGlvbn0gY2FsbGJhY2sgIGBjYWxsYmFjayh2YWx1ZSwga2V5KWBcclxuICAgICAqIEBwYXJhbSAgICAgICAge09iamVjdH0gICBjb250ZXh0ICAgIFRoZSBgdGhpc2AgYmluZGluZyBmb3IgY2FsbGJhY2tzLCBkZWZhdWx0IG51bGxcclxuICAgICAqL1xyXG4gICAgZnVuY3Rpb24gZm9yRWFjaChjYWxsYmFjaywgY29udGV4dCl7XHJcbiAgICAgIHZhciBkYXRhID0gdW53cmFwKHRoaXMpO1xyXG4gICAgICBjb250ZXh0ID0gY29udGV4dCA9PSBudWxsID8gZ2xvYmFsIDogY29udGV4dDtcclxuICAgICAgZm9yICh2YXIgaT0wOyBpIDwgMzsgaSsrKVxyXG4gICAgICAgIGZvciAodmFyIGtleSBpbiBkYXRhW2ldKVxyXG4gICAgICAgICAgY2FsbChjYWxsYmFjaywgY29udGV4dCwgZGF0YVtpXVtrZXldLCB1bmNvZXJjZShpLCBrZXkpLCB0aGlzKTtcclxuICAgIH1cclxuXHJcbiAgICBkZWxldGVfID0gZml4RGVsZXRlKGRlbGV0ZV8sIFsndmFsaWRhdGUnLCAndW53cmFwJywgJ2NvZXJjZSddLCBbdmFsaWRhdGUsIHVud3JhcCwgY29lcmNlXSk7XHJcbiAgICByZXR1cm4gW0hhc2hNYXAsIGdldCwgc2V0LCBoYXMsIGRlbGV0ZV8sIHNpemUsIGZvckVhY2hdO1xyXG4gIH0pO1xyXG5cclxuXHJcbiAgLy8gIyMjIyMjIyMjIyNcclxuICAvLyAjIyMgTWFwICMjI1xyXG4gIC8vICMjIyMjIyMjIyMjXHJcblxyXG4gIC8vIGlmIGEgZnVsbHkgaW1wbGVtZW50ZWQgTWFwIGV4aXN0cyB0aGVuIHVzZSBpdFxyXG4gIGlmICgnTWFwJyBpbiBnbG9iYWwgJiYgJ2ZvckVhY2gnIGluIGdsb2JhbC5NYXAucHJvdG90eXBlKSB7XHJcbiAgICBNID0gZXhwb3J0cy5NYXAgPSBnbG9iYWwuTWFwO1xyXG4gIH0gZWxzZSB7XHJcbiAgICBNID0gZXhwb3J0ZXIoJ01hcCcsIGZ1bmN0aW9uKHdyYXAsIHVud3JhcCl7XHJcbiAgICAgIC8vIGF0dGVtcHQgdG8gdXNlIGFuIGV4aXN0aW5nIHBhcnRpYWxseSBpbXBsZW1lbnRlZCBNYXBcclxuICAgICAgdmFyIEJ1aWx0aW5NYXAgPSBnbG9iYWwuTWFwLFxyXG4gICAgICAgICAgcHJvdG90eXBlID0gTWFwW3Byb3RvdHlwZV9dLFxyXG4gICAgICAgICAgd20gPSBXTVtwcm90b3R5cGVfXSxcclxuICAgICAgICAgIGhtID0gKEJ1aWx0aW5NYXAgfHwgSE0pW3Byb3RvdHlwZV9dLFxyXG4gICAgICAgICAgbWdldCAgICA9IFtjYWxsYmluZChobS5nZXQpLCBjYWxsYmluZCh3bS5nZXQpXSxcclxuICAgICAgICAgIG1zZXQgICAgPSBbY2FsbGJpbmQoaG0uc2V0KSwgY2FsbGJpbmQod20uc2V0KV0sXHJcbiAgICAgICAgICBtaGFzICAgID0gW2NhbGxiaW5kKGhtLmhhcyksIGNhbGxiaW5kKHdtLmhhcyldLFxyXG4gICAgICAgICAgbWRlbGV0ZSA9IFtjYWxsYmluZChobVsnZGVsZXRlJ10pLCBjYWxsYmluZCh3bVsnZGVsZXRlJ10pXTtcclxuXHJcbiAgICAgIHZhciB0eXBlID0gQnVpbHRpbk1hcFxyXG4gICAgICAgID8gZnVuY3Rpb24oKXsgcmV0dXJuIDAgfVxyXG4gICAgICAgIDogZnVuY3Rpb24obyl7IHJldHVybiArKHR5cGVvZiBvID09PSBvYmplY3RfID8gbyAhPT0gbnVsbCA6IHR5cGVvZiBvID09PSBmdW5jdGlvbl8pIH1cclxuXHJcbiAgICAgIC8vIGlmIHdlIGhhdmUgYSBidWlsdGluIE1hcCB3ZSBjYW4gbGV0IGl0IGRvIG1vc3Qgb2YgdGhlIGhlYXZ5IGxpZnRpbmdcclxuICAgICAgdmFyIGluaXQgPSBCdWlsdGluTWFwXHJcbiAgICAgICAgPyBmdW5jdGlvbigpeyByZXR1cm4geyAwOiBuZXcgQnVpbHRpbk1hcCB9IH1cclxuICAgICAgICA6IGZ1bmN0aW9uKCl7IHJldHVybiB7IDA6IG5ldyBITSwgMTogbmV3IFdNIH0gfTtcclxuXHJcbiAgICAgIC8qKlxyXG4gICAgICAgKiBAY2xhc3MgICAgICAgICBNYXBcclxuICAgICAgICogQGRlc2NyaXB0aW9uICAgQ29sbGVjdGlvbiB0aGF0IGFsbG93cyBhbnkga2luZCBvZiB2YWx1ZSB0byBiZSBhIGtleS5cclxuICAgICAgICogQHBhcmFtICAgICAgICAge0l0ZXJhYmxlfSBbaXRlcmFibGVdICBBbiBpdGVtIHRvIHBvcHVsYXRlIHRoZSBjb2xsZWN0aW9uIHdpdGguXHJcbiAgICAgICAqL1xyXG4gICAgICBmdW5jdGlvbiBNYXAoaXRlcmFibGUpe1xyXG4gICAgICAgIGlmICh0aGlzID09PSBnbG9iYWwgfHwgdGhpcyA9PSBudWxsIHx8IHRoaXMgPT09IHByb3RvdHlwZSlcclxuICAgICAgICAgIHJldHVybiBuZXcgTWFwKGl0ZXJhYmxlKTtcclxuXHJcbiAgICAgICAgdmFyIGRhdGEgPSBpbml0KCk7XHJcbiAgICAgICAgZGF0YS5rZXlzID0gW107XHJcbiAgICAgICAgZGF0YS52YWx1ZXMgPSBbXTtcclxuICAgICAgICB3cmFwKHRoaXMsIGRhdGEpO1xyXG5cclxuICAgICAgICB2YXIgc2VsZiA9IHRoaXM7XHJcbiAgICAgICAgaXRlcmFibGUgJiYgaW5pdGlhbGl6ZShpdGVyYWJsZSwgZnVuY3Rpb24odmFsdWUsIGtleSl7XHJcbiAgICAgICAgICBjYWxsKHNldCwgc2VsZiwgdmFsdWUsIGtleSk7XHJcbiAgICAgICAgfSk7XHJcbiAgICAgIH1cclxuICAgICAgLyoqXHJcbiAgICAgICAqIEBtZXRob2QgICAgICAgPGdldD5cclxuICAgICAgICogQGRlc2NyaXB0aW9uICBSZXRyaWV2ZSB0aGUgdmFsdWUgaW4gdGhlIGNvbGxlY3Rpb24gdGhhdCBtYXRjaGVzIGtleVxyXG4gICAgICAgKiBAcGFyYW0gICAgICAgIHtBbnl9IGtleVxyXG4gICAgICAgKiBAcmV0dXJuICAgICAgIHtBbnl9XHJcbiAgICAgICAqL1xyXG4gICAgICBmdW5jdGlvbiBnZXQoa2V5KXtcclxuICAgICAgICB2YXIgZGF0YSA9IHVud3JhcCh0aGlzKSxcclxuICAgICAgICAgICAgdCA9IHR5cGUoa2V5KTtcclxuICAgICAgICByZXR1cm4gZGF0YS52YWx1ZXNbbWdldFt0XShkYXRhW3RdLCBrZXkpXTtcclxuICAgICAgfVxyXG4gICAgICAvKipcclxuICAgICAgICogQG1ldGhvZCAgICAgICA8c2V0PlxyXG4gICAgICAgKiBAZGVzY3JpcHRpb24gIEFkZCBvciB1cGRhdGUgYSBwYWlyIGluIHRoZSBjb2xsZWN0aW9uLiBFbmZvcmNlcyB1bmlxdWVuZXNzIGJ5IG92ZXJ3cml0aW5nLlxyXG4gICAgICAgKiBAcGFyYW0gICAgICAgIHtBbnl9IGtleVxyXG4gICAgICAgKiBAcGFyYW0gICAgICAgIHtBbnl9IHZhbFxyXG4gICAgICAgKiovXHJcbiAgICAgIGZ1bmN0aW9uIHNldChrZXksIHZhbHVlKXtcclxuICAgICAgICB2YXIgZGF0YSA9IHVud3JhcCh0aGlzKSxcclxuICAgICAgICAgICAgdCA9IHR5cGUoa2V5KSxcclxuICAgICAgICAgICAgaW5kZXggPSBtZ2V0W3RdKGRhdGFbdF0sIGtleSk7XHJcblxyXG4gICAgICAgIGlmIChpbmRleCA9PT0gdW5kZWZpbmVkKSB7XHJcbiAgICAgICAgICBtc2V0W3RdKGRhdGFbdF0sIGtleSwgZGF0YS5rZXlzLmxlbmd0aCk7XHJcbiAgICAgICAgICBwdXNoKGRhdGEua2V5cywga2V5KTtcclxuICAgICAgICAgIHB1c2goZGF0YS52YWx1ZXMsIHZhbHVlKTtcclxuICAgICAgICB9IGVsc2Uge1xyXG4gICAgICAgICAgZGF0YS5rZXlzW2luZGV4XSA9IGtleTtcclxuICAgICAgICAgIGRhdGEudmFsdWVzW2luZGV4XSA9IHZhbHVlO1xyXG4gICAgICAgIH1cclxuICAgICAgfVxyXG4gICAgICAvKipcclxuICAgICAgICogQG1ldGhvZCAgICAgICA8aGFzPlxyXG4gICAgICAgKiBAZGVzY3JpcHRpb24gIENoZWNrIGlmIGtleSBleGlzdHMgaW4gdGhlIGNvbGxlY3Rpb24uXHJcbiAgICAgICAqIEBwYXJhbSAgICAgICAge0FueX0ga2V5XHJcbiAgICAgICAqIEByZXR1cm4gICAgICAge0Jvb2xlYW59IGlzIGluIGNvbGxlY3Rpb25cclxuICAgICAgICoqL1xyXG4gICAgICBmdW5jdGlvbiBoYXMoa2V5KXtcclxuICAgICAgICB2YXIgdCA9IHR5cGUoa2V5KTtcclxuICAgICAgICByZXR1cm4gbWhhc1t0XSh1bndyYXAodGhpcylbdF0sIGtleSk7XHJcbiAgICAgIH1cclxuICAgICAgLyoqXHJcbiAgICAgICAqIEBtZXRob2QgICAgICAgPGRlbGV0ZT5cclxuICAgICAgICogQGRlc2NyaXB0aW9uICBSZW1vdmUga2V5IGFuZCBtYXRjaGluZyB2YWx1ZSBpZiBmb3VuZFxyXG4gICAgICAgKiBAcGFyYW0gICAgICAgIHtBbnl9IGtleVxyXG4gICAgICAgKiBAcmV0dXJuICAgICAgIHtCb29sZWFufSB0cnVlIGlmIGl0ZW0gd2FzIGluIGNvbGxlY3Rpb25cclxuICAgICAgICovXHJcbiAgICAgIGZ1bmN0aW9uIGRlbGV0ZV8oa2V5KXtcclxuICAgICAgICB2YXIgZGF0YSA9IHVud3JhcCh0aGlzKSxcclxuICAgICAgICAgICAgdCA9IHR5cGUoa2V5KSxcclxuICAgICAgICAgICAgaW5kZXggPSBtZ2V0W3RdKGRhdGFbdF0sIGtleSk7XHJcblxyXG4gICAgICAgIGlmIChpbmRleCA9PT0gdW5kZWZpbmVkKVxyXG4gICAgICAgICAgcmV0dXJuIGZhbHNlO1xyXG5cclxuICAgICAgICBtZGVsZXRlW3RdKGRhdGFbdF0sIGtleSk7XHJcbiAgICAgICAgc3BsaWNlKGRhdGEua2V5cywgaW5kZXgsIDEpO1xyXG4gICAgICAgIHNwbGljZShkYXRhLnZhbHVlcywgaW5kZXgsIDEpO1xyXG4gICAgICAgIHJldHVybiB0cnVlO1xyXG4gICAgICB9XHJcbiAgICAgIC8qKlxyXG4gICAgICAgKiBAbWV0aG9kICAgICAgIDxzaXplPlxyXG4gICAgICAgKiBAZGVzY3JpcHRpb24gIFJldHJpZXZlIHRoZSBhbW91bnQgb2YgaXRlbXMgaW4gdGhlIGNvbGxlY3Rpb25cclxuICAgICAgICogQHJldHVybiAgICAgICB7TnVtYmVyfVxyXG4gICAgICAgKi9cclxuICAgICAgZnVuY3Rpb24gc2l6ZSgpe1xyXG4gICAgICAgIHJldHVybiB1bndyYXAodGhpcykua2V5cy5sZW5ndGg7XHJcbiAgICAgIH1cclxuICAgICAgLyoqXHJcbiAgICAgICAqIEBtZXRob2QgICAgICAgPGZvckVhY2g+XHJcbiAgICAgICAqIEBkZXNjcmlwdGlvbiAgTG9vcCB0aHJvdWdoIHRoZSBjb2xsZWN0aW9uIHJhaXNpbmcgY2FsbGJhY2sgZm9yIGVhY2hcclxuICAgICAgICogQHBhcmFtICAgICAgICB7RnVuY3Rpb259IGNhbGxiYWNrICBgY2FsbGJhY2sodmFsdWUsIGtleSlgXHJcbiAgICAgICAqIEBwYXJhbSAgICAgICAge09iamVjdH0gICBjb250ZXh0ICAgIFRoZSBgdGhpc2AgYmluZGluZyBmb3IgY2FsbGJhY2tzLCBkZWZhdWx0IG51bGxcclxuICAgICAgICovXHJcbiAgICAgIGZ1bmN0aW9uIGZvckVhY2goY2FsbGJhY2ssIGNvbnRleHQpe1xyXG4gICAgICAgIHZhciBkYXRhID0gdW53cmFwKHRoaXMpLFxyXG4gICAgICAgICAgICBrZXlzID0gZGF0YS5rZXlzLFxyXG4gICAgICAgICAgICB2YWx1ZXMgPSBkYXRhLnZhbHVlcztcclxuXHJcbiAgICAgICAgY29udGV4dCA9IGNvbnRleHQgPT0gbnVsbCA/IGdsb2JhbCA6IGNvbnRleHQ7XHJcblxyXG4gICAgICAgIGZvciAodmFyIGk9MCwgbGVuPWtleXMubGVuZ3RoOyBpIDwgbGVuOyBpKyspXHJcbiAgICAgICAgICBjYWxsKGNhbGxiYWNrLCBjb250ZXh0LCB2YWx1ZXNbaV0sIGtleXNbaV0sIHRoaXMpO1xyXG4gICAgICB9XHJcblxyXG4gICAgICBkZWxldGVfID0gZml4RGVsZXRlKGRlbGV0ZV8sXHJcbiAgICAgICAgWyd0eXBlJywgJ3Vud3JhcCcsICdjYWxsJywgJ3NwbGljZSddLFxyXG4gICAgICAgIFt0eXBlLCB1bndyYXAsIGNhbGwsIHNwbGljZV1cclxuICAgICAgKTtcclxuICAgICAgcmV0dXJuIFtNYXAsIGdldCwgc2V0LCBoYXMsIGRlbGV0ZV8sIHNpemUsIGZvckVhY2hdO1xyXG4gICAgfSk7XHJcbiAgfVxyXG5cclxuXHJcbiAgLy8gIyMjIyMjIyMjIyNcclxuICAvLyAjIyMgU2V0ICMjI1xyXG4gIC8vICMjIyMjIyMjIyMjXHJcblxyXG4gIGV4cG9ydGVyKCdTZXQnLCBmdW5jdGlvbih3cmFwLCB1bndyYXApe1xyXG4gICAgdmFyIHByb3RvdHlwZSA9IFNldFtwcm90b3R5cGVfXSxcclxuICAgICAgICBtID0gTVtwcm90b3R5cGVfXSxcclxuICAgICAgICBtc2l6ZSA9IGNhbGxiaW5kKG0uc2l6ZSksXHJcbiAgICAgICAgbWZvckVhY2ggPSBjYWxsYmluZChtLmZvckVhY2gpLFxyXG4gICAgICAgIG1nZXQgPSBjYWxsYmluZChtLmdldCksXHJcbiAgICAgICAgbXNldCA9IGNhbGxiaW5kKG0uc2V0KSxcclxuICAgICAgICBtaGFzID0gY2FsbGJpbmQobS5oYXMpLFxyXG4gICAgICAgIG1kZWxldGUgPSBjYWxsYmluZChtWydkZWxldGUnXSk7XHJcblxyXG4gICAgLyoqXHJcbiAgICAgKiBAY2xhc3MgICAgICAgIFNldFxyXG4gICAgICogQGRlc2NyaXB0aW9uICBDb2xsZWN0aW9uIG9mIHZhbHVlcyB0aGF0IGVuZm9yY2VzIHVuaXF1ZW5lc3MuXHJcbiAgICAgKiBAcGFyYW0gICAgICAgIHtJdGVyYWJsZX0gW2l0ZXJhYmxlXSAgQW4gaXRlbSB0byBwb3B1bGF0ZSB0aGUgY29sbGVjdGlvbiB3aXRoLlxyXG4gICAgICoqL1xyXG4gICAgZnVuY3Rpb24gU2V0KGl0ZXJhYmxlKXtcclxuICAgICAgaWYgKHRoaXMgPT09IGdsb2JhbCB8fCB0aGlzID09IG51bGwgfHwgdGhpcyA9PT0gcHJvdG90eXBlKVxyXG4gICAgICAgIHJldHVybiBuZXcgU2V0KGl0ZXJhYmxlKTtcclxuXHJcbiAgICAgIHdyYXAodGhpcywgbmV3IE0pO1xyXG5cclxuICAgICAgdmFyIHNlbGYgPSB0aGlzO1xyXG4gICAgICBpdGVyYWJsZSAmJiBpbml0aWFsaXplKGl0ZXJhYmxlLCBmdW5jdGlvbih2YWx1ZSwga2V5KXtcclxuICAgICAgICBjYWxsKGFkZCwgc2VsZiwga2V5KTtcclxuICAgICAgfSk7XHJcbiAgICB9XHJcbiAgICAvKipcclxuICAgICAqIEBtZXRob2QgICAgICAgPGFkZD5cclxuICAgICAqIEBkZXNjcmlwdGlvbiAgSW5zZXJ0IHZhbHVlIGlmIG5vdCBmb3VuZCwgZW5mb3JjaW5nIHVuaXF1ZW5lc3MuXHJcbiAgICAgKiBAcGFyYW0gICAgICAgIHtBbnl9IHZhbFxyXG4gICAgICovXHJcbiAgICBmdW5jdGlvbiBhZGQoa2V5KXtcclxuICAgICAgbXNldCh1bndyYXAodGhpcyksIGtleSwga2V5KTtcclxuICAgIH1cclxuICAgIC8qKlxyXG4gICAgICogQG1ldGhvZCAgICAgICA8aGFzPlxyXG4gICAgICogQGRlc2NyaXB0aW9uICBDaGVjayBpZiBrZXkgZXhpc3RzIGluIHRoZSBjb2xsZWN0aW9uLlxyXG4gICAgICogQHBhcmFtICAgICAgICB7QW55fSBrZXlcclxuICAgICAqIEByZXR1cm4gICAgICAge0Jvb2xlYW59IGlzIGluIGNvbGxlY3Rpb25cclxuICAgICAqKi9cclxuICAgIGZ1bmN0aW9uIGhhcyhrZXkpe1xyXG4gICAgICByZXR1cm4gbWhhcyh1bndyYXAodGhpcyksIGtleSk7XHJcbiAgICB9XHJcbiAgICAvKipcclxuICAgICAqIEBtZXRob2QgICAgICAgPGRlbGV0ZT5cclxuICAgICAqIEBkZXNjcmlwdGlvbiAgUmVtb3ZlIGtleSBhbmQgbWF0Y2hpbmcgdmFsdWUgaWYgZm91bmRcclxuICAgICAqIEBwYXJhbSAgICAgICAge0FueX0ga2V5XHJcbiAgICAgKiBAcmV0dXJuICAgICAgIHtCb29sZWFufSB0cnVlIGlmIGl0ZW0gd2FzIGluIGNvbGxlY3Rpb25cclxuICAgICAqL1xyXG4gICAgZnVuY3Rpb24gZGVsZXRlXyhrZXkpe1xyXG4gICAgICByZXR1cm4gbWRlbGV0ZSh1bndyYXAodGhpcyksIGtleSk7XHJcbiAgICB9XHJcbiAgICAvKipcclxuICAgICAqIEBtZXRob2QgICAgICAgPHNpemU+XHJcbiAgICAgKiBAZGVzY3JpcHRpb24gIFJldHJpZXZlIHRoZSBhbW91bnQgb2YgaXRlbXMgaW4gdGhlIGNvbGxlY3Rpb25cclxuICAgICAqIEByZXR1cm4gICAgICAge051bWJlcn1cclxuICAgICAqL1xyXG4gICAgZnVuY3Rpb24gc2l6ZSgpe1xyXG4gICAgICByZXR1cm4gbXNpemUodW53cmFwKHRoaXMpKTtcclxuICAgIH1cclxuICAgIC8qKlxyXG4gICAgICogQG1ldGhvZCAgICAgICA8Zm9yRWFjaD5cclxuICAgICAqIEBkZXNjcmlwdGlvbiAgTG9vcCB0aHJvdWdoIHRoZSBjb2xsZWN0aW9uIHJhaXNpbmcgY2FsbGJhY2sgZm9yIGVhY2guIEluZGV4IGlzIHNpbXBseSB0aGUgY291bnRlciBmb3IgdGhlIGN1cnJlbnQgaXRlcmF0aW9uLlxyXG4gICAgICogQHBhcmFtICAgICAgICB7RnVuY3Rpb259IGNhbGxiYWNrICBgY2FsbGJhY2sodmFsdWUsIGluZGV4KWBcclxuICAgICAqIEBwYXJhbSAgICAgICAge09iamVjdH0gICBjb250ZXh0ICAgIFRoZSBgdGhpc2AgYmluZGluZyBmb3IgY2FsbGJhY2tzLCBkZWZhdWx0IG51bGxcclxuICAgICAqL1xyXG4gICAgZnVuY3Rpb24gZm9yRWFjaChjYWxsYmFjaywgY29udGV4dCl7XHJcbiAgICAgIHZhciBpbmRleCA9IDAsXHJcbiAgICAgICAgICBzZWxmID0gdGhpcztcclxuICAgICAgbWZvckVhY2godW53cmFwKHRoaXMsIGZ1bmN0aW9uKGtleSl7XHJcbiAgICAgICAgY2FsbChjYWxsYmFjaywgdGhpcywga2V5LCBpbmRleCsrLCBzZWxmKTtcclxuICAgICAgfSwgY29udGV4dCkpO1xyXG4gICAgfVxyXG5cclxuICAgIGRlbGV0ZV8gPSBmaXhEZWxldGUoZGVsZXRlXywgWydtZGVsZXRlJywgJ3Vud3JhcCddLCBbbWRlbGV0ZSwgdW53cmFwXSk7XHJcbiAgICByZXR1cm4gW1NldCwgYWRkLCBoYXMsIGRlbGV0ZV8sIHNpemUsIGZvckVhY2hdO1xyXG4gIH0pO1xyXG59KCdzdHJpbmcnLCAnb2JqZWN0JywgJ2Z1bmN0aW9uJywgJ3Byb3RvdHlwZScsICd0b1N0cmluZycsXHJcbiAgQXJyYXksIE9iamVjdCwgRnVuY3Rpb24sIEZ1bmN0aW9uLnByb3RvdHlwZSwgKDAsIGV2YWwpKCd0aGlzJyksXHJcbiAgdHlwZW9mIGV4cG9ydHMgPT09ICd1bmRlZmluZWQnID8gdGhpcyA6IGV4cG9ydHMsIHt9KTtcclxuIiwiXG5tb2R1bGUuZXhwb3J0cyA9IHBhcnNlXG5cbi8qKlxuICogZXhwZWN0ZWQgYXJndW1lbnQgbGVuZ3Roc1xuICogQHR5cGUge09iamVjdH1cbiAqL1xuXG52YXIgbGVuZ3RoID0ge2E6IDcsIGM6IDYsIGg6IDEsIGw6IDIsIG06IDIsIHE6IDQsIHM6IDQsIHQ6IDIsIHY6IDEsIHo6IDB9XG5cbi8qKlxuICogc2VnbWVudCBwYXR0ZXJuXG4gKiBAdHlwZSB7UmVnRXhwfVxuICovXG5cbnZhciBzZWdtZW50ID0gLyhbYXN0dnpxbWhsY10pKFteYXN0dnpxbWhsY10qKS9pZ1xuXG4vKipcbiAqIHBhcnNlIGFuIHN2ZyBwYXRoIGRhdGEgc3RyaW5nLiBHZW5lcmF0ZXMgYW4gQXJyYXlcbiAqIG9mIGNvbW1hbmRzIHdoZXJlIGVhY2ggY29tbWFuZCBpcyBhbiBBcnJheSBvZiB0aGVcbiAqIGZvcm0gYFtjb21tYW5kLCBhcmcxLCBhcmcyLCAuLi5dYFxuICpcbiAqIEBwYXJhbSB7U3RyaW5nfSBwYXRoXG4gKiBAcmV0dXJuIHtBcnJheX1cbiAqL1xuXG5mdW5jdGlvbiBwYXJzZShwYXRoKSB7XG5cdHZhciBkYXRhID0gW11cblx0cGF0aC5yZXBsYWNlKHNlZ21lbnQsIGZ1bmN0aW9uKF8sIGNvbW1hbmQsIGFyZ3Mpe1xuXHRcdHZhciB0eXBlID0gY29tbWFuZC50b0xvd2VyQ2FzZSgpXG5cdFx0YXJncyA9IHBhcnNlVmFsdWVzKGFyZ3MpXG5cblx0XHQvLyBvdmVybG9hZGVkIG1vdmVUb1xuXHRcdGlmICh0eXBlID09ICdtJyAmJiBhcmdzLmxlbmd0aCA+IDIpIHtcblx0XHRcdGRhdGEucHVzaChbY29tbWFuZF0uY29uY2F0KGFyZ3Muc3BsaWNlKDAsIDIpKSlcblx0XHRcdHR5cGUgPSAnbCdcblx0XHRcdGNvbW1hbmQgPSBjb21tYW5kID09ICdtJyA/ICdsJyA6ICdMJ1xuXHRcdH1cblxuXHRcdHdoaWxlICh0cnVlKSB7XG5cdFx0XHRpZiAoYXJncy5sZW5ndGggPT0gbGVuZ3RoW3R5cGVdKSB7XG5cdFx0XHRcdGFyZ3MudW5zaGlmdChjb21tYW5kKVxuXHRcdFx0XHRyZXR1cm4gZGF0YS5wdXNoKGFyZ3MpXG5cdFx0XHR9XG5cdFx0XHRpZiAoYXJncy5sZW5ndGggPCBsZW5ndGhbdHlwZV0pIHRocm93IG5ldyBFcnJvcignbWFsZm9ybWVkIHBhdGggZGF0YScpXG5cdFx0XHRkYXRhLnB1c2goW2NvbW1hbmRdLmNvbmNhdChhcmdzLnNwbGljZSgwLCBsZW5ndGhbdHlwZV0pKSlcblx0XHR9XG5cdH0pXG5cdHJldHVybiBkYXRhXG59XG5cbmZ1bmN0aW9uIHBhcnNlVmFsdWVzKGFyZ3Mpe1xuXHRhcmdzID0gYXJncy5tYXRjaCgvLT9bLjAtOV0rKD86ZVstK10/XFxkKyk/L2lnKVxuXHRyZXR1cm4gYXJncyA/IGFyZ3MubWFwKE51bWJlcikgOiBbXVxufVxuIiwiXCJ1c2Ugc3RyaWN0XCJcblxudmFyIHNoZWxsID0gcmVxdWlyZShcImdhbWUtc2hlbGxcIikoKTtcbnZhciBNYXAgPSByZXF1aXJlKFwiaGFybW9ueS1jb2xsZWN0aW9uc1wiKS5NYXA7XG52YXIgU2V0ID0gcmVxdWlyZShcImhhcm1vbnktY29sbGVjdGlvbnNcIikuU2V0O1xudmFyIGR0ID0gcmVxdWlyZShcImRlbGF1bmF5LXRyaWFuZ3VsYXRlXCIpO1xudmFyIHBhcnNlID0gcmVxdWlyZSgncGFyc2Utc3ZnLXBhdGgnKVxuXG5cbnZhciBuYlJhbmRvbVBvaW50cyA9IDEwMDtcbnZhciBuYkFudHMgPSAyMDA7XG52YXIgdGV4dE1lc2ggPSBmYWxzZTtcbnZhciBuYkNpdHkgPSAxMDtcblxudmFyIHNxcnQgPSBNYXRoLnNxcnQ7XG52YXIgcG93ID0gTWF0aC5wb3c7XG52YXIgZmxvb3IgPSBNYXRoLmZsb29yO1xudmFyIHJhbmRvbSA9IE1hdGgucmFuZG9tO1xuXG52YXIgc3ZnU3RyaW5nID0gXCJtIDEyNDYuMzg2NCw0MjIuNjcwNDYgYyAtMzIuNjU2NiwtMy4wNTk3NyAtODAuOTIwMiwtMS42MDg5MiAtNjQsLTQ5LjIzNTcyIC0xLjQ3NTksLTE3LjI2OTk5IC0yLjkzODksLTM0LjUyNjk3IDE5LjYxNzksLTI3Ljc2NDI4IDI4LjIwNTEsLTcuNzE3ODIgNS42Mjk5LDM2LjUzNDI4IDI5LjM4MjEsNDQgMTUuNjQwMSwxMS4wODc4OSAzNi44MjY4LDE1LjE2Nzg4IDU1Ljc2MDEsMTEuNDM1OTEgMTguOTMzMywtMy43MzE5OCAzNS42MTMyLC0xNS4yNzU5MyA0Mi4yMzk5LC0zNS40MzU5MSAxMS4xMTQsLTMwLjUzNDQ3IC01LjQ1NSwtNjMuNjY5OTcgLTM0Ljc3OTIsLTc2IC0yMS4yMzc5LC0xMi4xODg3NCAtNDUuMTk2MywtMjQuMzE2NDcgLTY1LjUxODYsLTM5LjgzMzE1IC0yMC4zMjIzLC0xNS41MTY2OCAtMzcuMDA4NSwtMzQuNDIyMzEgLTQzLjcwMjIsLTYwLjE2Njg1IC01LjI3NzcsLTIzLjgyNjUxIDEuMjk1NywtNDcuNDgxNTEgMTUuMDI2OSwtNjYuMTAyMzggMTMuNzMxMiwtMTguNjIwODggMzQuNjIwMSwtMzIuMjA3NjIzIDU3Ljk3MzEsLTM1Ljg5NzYxOSAxOS41NzE0LC00LjQwODA2NyAzOS41MTk5LC00LjQ1NTIxNyA1OS4zMTMyLC0yLjA1NTU4NCAxOS43OTMzLDIuMzk5NjM0IDM5LjQzMTMsNy4yNDYwNTEgNTguMzgxNywxMi42MjUxMiAxLjUwOTYsMTkuMTczMDIzIDUuNzMyOSw0OS42NzAxMzMgLTAuMzEyOCw2NS40MzA0NjMgLTE2LjY1NywtMC4yNjkxNSAtMzkuMTg0Miw1Ljg4NDE0IC0zMS42NTQ4LC0xOSAtMS42NDgyLC0yMy44MjA1MyAtMjEuMzQyNywtMzYuODI3MDEgLTQyLjQ5MzgsLTM4LjUwMDEyIC0yMS4xNTExLC0xLjY3MzEyIC00My43NTg4LDcuOTg3MTUgLTUxLjIzMzUsMjkuNTAwMTIgLTE0LjE4MjYsMjcuODQ2MDMgNC4xOTE0LDU3LjkyMDIgMzAsNzAgMTkuOTY0OSwxMy40OTQyMyA0My41NTI2LDI0Ljc2NDM4IDY0LjM1NTIsMzguNzc1ODkgMjAuODAyNywxNC4wMTE1MiAzOC44MjAyLDMwLjc2NDQgNDcuNjQ0OCw1NS4yMjQxMSA4Ljg1NzcsMjQuNzM2NjcgNC43MzIxLDUwLjkwMDgyIC03Ljc5NTEsNzIuNDIzNzIgLTEyLjUyNzIsMjEuNTIyOSAtMzMuNDU2LDM4LjQwNDU2IC01OC4yMDQ5LDQ0LjU3NjI4IC0yNS41NTE2LDguODEzNTkgLTUzLjQyODUsNy42MTM5MSAtODAsNiB6IG0gLTUwMi43NTAwNCwtNC40OTU4MyBjIC0xNC40OTYzNCwtMjAuMTY3OTggLTI5LjA1MTc3LC00MC4zMTk4OSAtNDMuNjY0OTQsLTYwLjQ0NDE2IC0xNC42MTMxNiwtMjAuMTI0MjggLTI5LjI4NDA0LC00MC4yMjA5MSAtNDQuMDExMjYsLTYwLjI3ODM0IC0xNC43MjcyMywtMjAuMDU3NDQgLTI5LjUxMDc5LC00MC4wNzU2NyAtNDQuMzQ5MzEsLTYwLjA0MzEzIC0xNC44Mzg1MSwtMTkuOTY3NDcgLTI5LjczMTk5LC0zOS44ODQxNyAtNDQuNjc5MDQsLTU5LjczODU0IDAuMDg1LDY3Ljk0MTE0IC0yLjU5NDU0LDEzNi4yODYzNSAxLjQ1NDU1LDIwNCAyLjM3MzYxLDIxLjMzMjk3IDQxLjk1OTY4LDIuNjA4MDIgMzIsMjkuMjM1NzIgLTUuNjg0ODQsMTAuNzM3NzggLTMwLjk3OTE4LDEuODY2MDggLTQ0LjUwMzYsNC43NjQyOCAtMTkuMTY1NDcsMCAtMzguMzMwOTMsMCAtNTcuNDk2NCwwIC0xMi4xNzc0NSwtMzEuOTk5ODMgMzIuMzQxMjYsLTExLjA5MzUxIDM0LC0zOCAwLjkyMjcsLTE5Ljk1MTQ3IDEuNTA2MTgsLTM5LjkzMDg3IDEuODQ3NDcsLTU5LjkyNjM0IDAuMzQxMjksLTE5Ljk5NTQ3IDAuNDQwMzksLTQwLjAwNzAyIDAuMzk0MzUsLTYwLjAyMjgxIC0wLjA0NjEsLTIwLjAxNTc5IC0wLjIzNzI1LC00MC4wMzU4MiAtMC40NzY1NywtNjAuMDQ4MjQgLTAuMjM5MzEsLTIwLjAxMjQyIC0wLjUyNjc0LC00MC4wMTcyNCAtMC43NjUyNSwtNjAuMDAyNjEgMTAuMjAyOTksLTMyLjgxMzk2IC02My44ODAzOSwtMTguNjY5MzIgLTI1LjA2NTY2LC00NS43MzE3MzkgMjMuODQ3MjQsMS4zMTMzMTIgNTAuMzM2NjksLTIuOTc1NDE2IDczLjAwNzY4LDEuNTgyODQ2IDE0LjU2NzM5LDIwLjI0MjQ1MyAyOS4yMzg2MSw0MC40MTk2NDMgNDMuOTU1MzEsNjAuNTY4NjMzIDE0LjcxNjcxLDIwLjE0ODk5IDI5LjQ3ODg5LDQwLjI2OTc3IDQ0LjIyODIsNjAuMzk5NDMgMTQuNzQ5MzEsMjAuMTI5NjUgMjkuNDg1NzQsNDAuMjY4MTcgNDQuMTUwOTQsNjAuNDUyNjIgMTQuNjY1MTksMjAuMTg0NDUgMjkuMjU5MTYsNDAuNDE0ODMgNDMuNzIzNTMsNjAuNzI4MjEgMC44MzMzMiwtNzAuMDcwMDIgMi44MjE0NSwtMTQwLjc4NjkyIC0xLC0yMTAuNTQ1NDYgLTEwLjM5NjAxLC05LjcxMjQ2IC01MS41NDAxNCwtMjEuODU5NjUgLTI0LjE0NjQyLC0zMy40NTQ1MzkgMzAuOTAwMywwLjQ5MTY3IDYxLjkyNjk2LC0xLjI5Mzg0MSA5Mi43MTc4NSwxLjQyODU3IDEyLjQ0OTE4LDMzLjYzODMxOSAtNTEuNjYwMDMsMTIuMjIzNzM5IC0zNC4wMDY5Miw1NS40MTUxNzkgLTAuMjQxMjEsMjIuNTk1MzIgLTAuMzg4MjQsNDUuMTkxMTcgLTAuNDc1MzIsNjcuNzg3MzQgLTAuMDg3MSwyMi41OTYxOCAtMC4xMTQxOSw0NS4xOTI2NyAtMC4xMTU1OSw2Ny43ODkzIC0xMGUtNCwyMi41OTY2MiAwLjAyMjksNDUuMTkzMzYgMC4wMzg3LDY3Ljc5MDAzIDAuMDE1OCwyMi41OTY2NyAwLjAyMzEsNDUuMTkzMjYgLTAuMDEyMyw2Ny43ODk1OCAtMTUuNTM1NiwtMC41Mzg4MyAtMzEuNDM3MjIsMS4yMzA5OCAtNDYuNzUsLTEuNDk1ODMgeiBNIDEzMi4zODYzNyw0MDUuNjcwNDYgYyAtMi4zMDE5NSwtMTkuMzgwOTggMzUuMjA4NTMsLTYuMDE4ODUgMzUuOTk5OTksLTMyIDEwLjY5NDMsLTIzLjIwMzQxIDIxLjM1ODA2LC00Ni40MTgzMyAzMi4wMDg5OSwtNjkuNjM4MTggMTAuNjUwOTMsLTIzLjIxOTg1IDIxLjI4OTA1LC00Ni40NDQ2MSAzMS45MzIwNywtNjkuNjY3NyAxMC42NDMwMywtMjMuMjIzMDggMjEuMjkwOTYsLTQ2LjQ0NDQ5IDMxLjk2MTU0LC02OS42NTc2MSAxMC42NzA1OCwtMjMuMjEzMTMgMjEuMzYzOCwtNDYuNDE3OTcgMzIuMDk3NCwtNjkuNjA3OTM5IDE0LjAwNjMzLC0xNC43ODQ2NjQgMzEuMDYzMDEsLTMuNTg5NTA4IDM0LDE0LjU3MTQyOSA5LjQ3NTg5LDIwLjUwOTk4IDE4LjkwMjY3LDQxLjA0MDMyIDI4LjMyMTY4LDYxLjU3Mzg2IDkuNDE5MDEsMjAuNTMzNTUgMTguODMwMjcsNDEuMDcwMzEgMjguMjc1MTIsNjEuNTkzMTUgOS40NDQ4NSwyMC41MjI4NCAxOC45MjMzMSw0MS4wMzE3NyAyOC40NzY3Miw2MS41MDk2NCA5LjU1MzQxLDIwLjQ3Nzg4IDE5LjE4MTc5LDQwLjkyNDcgMjguOTI2NDgsNjEuMzIzMzUgOC45NDg0NywxOC41NzQ5NCAxNy45NzY4NCw0NC41OTk3NiA0NCw0MCA2Ljk1MjgzLDIzLjIxODMxIC0xMS4yMjYxOSwyMS4zMDQ3NyAtMjguNzk2NDUsMjAgLTE3LjUzMzkzLDAgLTM1LjA2Nzg1LDAgLTUyLjYwMTc4LDAgLTE3LjUzMzkyLDAgLTM1LjA2Nzg1LDAgLTUyLjYwMTc3LDAgLTEzLjYxMDUzLC0zMS41Njk1MyAzNC44NjQ3OCwtMTEuNDMwOTkgMzAsLTM3IC03Ljk2MzI2LC0yMC40NDkzNyAtMTYuOTAwNjcsLTQwLjg4MjE5IC0yNy4xNzE4OCwtNjAuMTQ4MSAtMjIuNjk4MzUsLTAuMzAxNTkgLTQ1Ljc0NzgzLC0xLjAwNDc4IC02OC43MDE0OSwtMS4xMzc0NiAtMjIuOTUzNjcsLTAuMTMyNjkgLTQ1LjgxMTUzLDAuMzA1MTMgLTY4LjEyNjYzLDIuMjg1NTYgLTEwLjE1NjE0LDIzLjQ4Mzk0IC0yNy42NTY1Miw0Ni42MjkyMSAtMjYuNTcxNDMsNzIuNTcxNDMgMjEuODc5ODMsLTkuNTk1MjkgMzcuNDQyMzcsMjMuODI3NjUgMTQuMDk5OTksMjMuNDI4NTcgLTI1LjE3NjE4LDAgLTUwLjM1MjM3LDAgLTc1LjUyODU1LDAgMCwtMy4zMzMzMyAwLC02LjY2NjY3IDAsLTEwIHogbSAyMTMuOTk5OTksLTExMiBjIC05Ljg1MTg5LC0xOS45MjkxMyAtMTkuMTI5NjgsLTQwLjEwODY0IC0yOC40NjQ5LC02MC4yNjAzNyAtOS4zMzUyMSwtMjAuMTUxNzMgLTE4LjcyNzg1LC00MC4yNzU2NyAtMjguODA5NCwtNjAuMDkzNjYgLTguMDQ3MSw0LjYyNTgyIC0xMy4xMjEwOCwyNi45NzYwNyAtMTkuNzI1NywzOC4zNTQwMyAtMTIuNjA4MjksMjcuMjE5NDQgLTI1LjI2NjYxLDU0LjQxODE2IC0zNyw4MiAxOC4xOTM0NiwxLjY1ODQ2IDM3LjYwMDA5LDIuNDg3NTcgNTcuMDA1MDUsMi40ODc1MSAxOS40MDQ5NSwtNmUtNSAzOC44MDgyMiwtMC44MjkyOSA1Ni45OTQ5NSwtMi40ODc1MSB6IG0gNTgyLDExMiBjIC0xLjIzMzcsLTIxLjg4MTg0IDUyLjc4NjQyLC0zLjk0MTAzIDM5LjQzNTMyLC00Mi4xNTYyNSAwLjI0MTk1LC0yMC44MTkyIDAuMzg5MjUsLTQxLjYzODk3IDAuNDc2MzUsLTYyLjQ1OTA5IDAuMDg3MSwtMjAuODIwMTIgMC4xMTM5OCwtNDEuNjQwNTggMC4xMTUxLC02Mi40NjExOCAxMGUtNCwtMjAuODIwNiAtMC4wMjM1LC00MS42NDEzNCAtMC4wMzk0LC02Mi40NjE5OSAtMC4wMTU5LC0yMC44MjA2NSAtMC4wMjMyLC00MS42NDEyMiAwLjAxMjcsLTYyLjQ2MTQ5IC0yNi4zNzgzMiwyLjI5ODI4IC01NS45MzIzMiwtNi4yMDM3MSAtNzkuODU3ODcsNi44MDU4OSAtMS41MDI1NiwyOC4wNDQ3MyAtNC42ODM0NCw0Ni42NDY1MSAtMzUuMjMxNDIsMzkuMTk0MTEgLTcuMTE2NDMsLTE1LjU0MTE2IC0wLjkwMTA0LC00Ni44MzU4NCAtMS40ODIxNCwtNjYuNTcxNDI5IDI1LjAyNDk1LC0wLjcwODEyIDUwLjA2MTcsLTEuMTE1Nzc5IDc1LjEwNTY1LC0xLjMzNzg3MyAyNS4wNDM5NSwtMC4yMjIwOTMgNTAuMDk1MSwtMC4yNTg2MjEgNzUuMTQ4ODUsLTAuMjI0NDc2IDI1LjA1MzgsMC4wMzQxNCA1MC4xMTAxLDAuMTM4OTYxIDc1LjE2NDUsMC4xOTk1NTUgMjUuMDU0NCwwLjA2MDU5IDUwLjEwNjcsMC4wNzY5NyA3NS4xNTI1LC0wLjA2NTc4IDAsMjIuNjY2NjczIDAsNDUuMzMzMzMzIDAsNjguMDAwMDAzIC0xMS41NzA4LC0zLjA5NTg0IC0zOC4xMTI4LDguNzc5NDUgLTMyLC0xMiA2LjM3NjMsLTQxLjUzODgzIC00MC4yOTQzLC0zNC4zNjI2MSAtNjcuNzA3MiwtMzQgLTIzLjYwMDIsLTcuNjk3ODQgLTE0Ljk0MDcsMTYuODY0MDkgLTE2LjI5MjksMzEuNDE0MzEgMC4xNzQ1LDE5Ljc1NzM5IDAuMTQ2NywzOS41MzAyNCAwLjA2Nyw1OS4zMDc2MyAtMC4wOCwxOS43NzczOCAtMC4yMTE2LDM5LjU1OTMgLTAuMjQ1Miw1OS4zMzQ4MSAtMC4wMzQsMTkuNzc1NTEgMC4wMzEsMzkuNTQ0NjIgMC4zNDM5LDU5LjI5NjM4IDAuMzEyOSwxOS43NTE3NiAwLjg3NDMsMzkuNDg2MTggMS44MzQzLDU5LjE5MjMyIDYuMzIwNCwyMi43Njc1NiA1MC4zOTgzLDAuMDU3NCAzOC4wMDAxLDMzLjQ1NDU1IC0yNC42NjY3LDAgLTQ5LjMzMzMsMCAtNzQsMCAtMjQuNjY2NzEsMCAtNDkuMzMzNCwwIC03NC4wMDAwNywwIDAsLTMuMzMzMzMgMCwtNi42NjY2NyAwLC0xMCB6XCI7XG5cbmZ1bmN0aW9uIHN2Z1RvUG9pbnRzKHN2Z1N0cmluZykge1xuICAgIHZhciBwb2ludHMgPSBbXTtcbiAgICB2YXIgZWRnZXMgPSBbXTtcblxuICAgIHZhciBYID0gMDtcbiAgICB2YXIgWSA9IDA7XG4gICAgdmFyIG5iUG9pbnRzID0gMDtcbiAgICB2YXIgcHJldlBvaW50O1xuICAgIHZhciBjb21tYW5kcyA9IHBhcnNlKHN2Z1N0cmluZylcbiAgICBmb3IgKHZhciBpPTA7IGk8Y29tbWFuZHMubGVuZ3RoOyBpKyspe1xuICAgICAgICB2YXIgY29tbWFuZCA9IGNvbW1hbmRzW2ldO1xuICAgICAgICBzd2l0Y2ggKGNvbW1hbmRbMF0pIHtcbiAgICAgICAgICAgIGNhc2UgXCJtXCI6XG4gICAgICAgICAgICAgICAgWCArPSBjb21tYW5kWzFdO1xuICAgICAgICAgICAgICAgIFkgKz0gY29tbWFuZFsyXTtcbiAgICAgICAgICAgICAgICBwcmV2UG9pbnQgPSB1bmRlZmluZWQ7XG4gICAgICAgICAgICAgICAgYnJlYWs7XG4gICAgICAgICAgICBjYXNlIFwiTVwiOlxuICAgICAgICAgICAgICAgIFggPSBjb21tYW5kWzFdO1xuICAgICAgICAgICAgICAgIFkgPSBjb21tYW5kWzJdO1xuICAgICAgICAgICAgICAgIHByZXZQb2ludCA9IHVuZGVmaW5lZDtcbiAgICAgICAgICAgICAgICBicmVhazsgIFxuICAgICAgICAgICAgY2FzZSBcImNcIjpcbiAgICAgICAgICAgICAgICBYICs9IGNvbW1hbmRbNV07XG4gICAgICAgICAgICAgICAgWSArPSBjb21tYW5kWzZdO1xuICAgICAgICAgICAgICAgIHBvaW50cy5wdXNoKHtpZDpuYlBvaW50cywgeDpYLCB5Oll9KTtcbiAgICAgICAgICAgICAgICBuYlBvaW50cysrO1xuICAgICAgICAgICAgICAgIGlmIChwcmV2UG9pbnQpIHtcbiAgICAgICAgICAgICAgICAgICAgZWRnZXMucHVzaChbcHJldlBvaW50LCBuYlBvaW50c10pO1xuICAgICAgICAgICAgICAgICAgICBwcmV2UG9pbnQgPSBuYlBvaW50cztcbiAgICAgICAgICAgICAgICB9XG4gICAgICAgICAgICAgICAgYnJlYWs7ICAgIFxuICAgICAgICB9XG4gICAgfVxuICAgIHJldHVybiB7cG9pbnRzIDogcG9pbnRzLCBlZGdlcyA6IGVkZ2VzfTtcbn1cblxuZnVuY3Rpb24gcmFuZ2Uoc3RhcnQsIGNvdW50KSB7XG4gICAgcmV0dXJuIEFycmF5LmFwcGx5KDAsIEFycmF5KGNvdW50KSkubWFwKGZ1bmN0aW9uIChlbGVtZW50LCBpbmRleCkgeyByZXR1cm4gaW5kZXggKyBzdGFydCB9KTtcbn1cblxuZnVuY3Rpb24gc2lnbih4KSB7IHJldHVybiB4ID8geCA8IDAgPyAtMSA6IDEgOiAwOyB9XG5cbi8vIGluaXRpYWxpemUgcG9pbnRzXG52YXIgcG9pbnRzID0gW107XG52YXIgY2l0eVNldDtcblxuaWYgKHRleHRNZXNoKXtcblxuICAgIHZhciBteVRleHQgPSBzdmdUb1BvaW50cyhzdmdTdHJpbmcpO1xuICAgIHBvaW50cyA9IG15VGV4dC5wb2ludHM7XG4gICAgY2l0eVNldCA9IG5ldyBTZXQocmFuZ2UoMCwgcG9pbnRzLmxlbmd0aCkpO1xuICAgIHZhciBzY2FsZSA9IDAuNVxuXG4gICAgLy8gc2NhbGUgcG9pbnRzIHRvIFswLDFdICsgc2NhbGVcbiAgICB2YXIgbWF4WCA9IE1hdGgubWF4LmFwcGx5KE1hdGgsIHBvaW50cy5tYXAoZnVuY3Rpb24ocCl7cmV0dXJuIHAueH0pKTtcbiAgICB2YXIgbWluWCA9IE1hdGgubWluLmFwcGx5KE1hdGgsIHBvaW50cy5tYXAoZnVuY3Rpb24ocCl7cmV0dXJuIHAueH0pKTtcbiAgICB2YXIgbWF4WSA9IE1hdGgubWF4LmFwcGx5KE1hdGgsIHBvaW50cy5tYXAoZnVuY3Rpb24ocCl7cmV0dXJuIHAueX0pKTtcbiAgICB2YXIgbWluWSA9IE1hdGgubWluLmFwcGx5KE1hdGgsIHBvaW50cy5tYXAoZnVuY3Rpb24ocCl7cmV0dXJuIHAueX0pKTtcbiAgICBwb2ludHMgPSBwb2ludHMubWFwKGZ1bmN0aW9uKHApe3JldHVybiB7aWQ6cC5pZCwgeDogMC40KihwLngtbWluWCkvKG1heFgtbWluWCkrMC4yNSwgeTogMC40KihwLnktbWluWSkvKG1heFktbWluWSkrMC4yNX19KVxuXG4gICAgLy9hZGQgcmFuZG9tIHBvaW50c1xuICAgIHZhciBuYlBvaW50cyA9IHBvaW50cy5sZW5ndGg7XG4gICAgZm9yKHZhciBpPTA7IGk8bmJSYW5kb21Qb2ludHM7ICsraSkge1xuICAgICAgICBwb2ludHMucHVzaCh7aWQgOiBuYlBvaW50cywgeDpyYW5kb20oKSwgeTpyYW5kb20oKX0pO1xuICAgICAgICBuYlBvaW50cysrO1xuICAgIH1cblxufSBlbHNlIHtcbiAgICAvL2FkZCByYW5kb20gcG9pbnRzXG4gICAgdmFyIG5iUG9pbnRzID0gMDtcbiAgICBmb3IodmFyIGk9MDsgaTxuYlJhbmRvbVBvaW50czsgKytpKSB7XG4gICAgICAgIHBvaW50cy5wdXNoKHtpZCA6IG5iUG9pbnRzLCB4OnJhbmRvbSgpLCB5OnJhbmRvbSgpfSk7XG4gICAgICAgIG5iUG9pbnRzKys7XG4gICAgfVxuICAgIGNpdHlTZXQgPSBuZXcgU2V0KHJhbmdlKDAsIG5iQ2l0eSkpO1xufVxuXG5cblxuXG4vLyB0cmlhbmd1bGF0ZVxudmFyIGNlbGxzID0gZHQocG9pbnRzLm1hcChmdW5jdGlvbihwKXtyZXR1cm4gW3AueCwgcC55XX0pKVxuXG4vLyBjcmVhdGUgZWRnZXNcbnZhciBuZXh0RWRnZXMgPSBuZXcgTWFwKCk7XG52YXIgZWRnZXMgPSBbXTtcbnZhciBwZXJtdXRhdGlvbnMgPSBbWzAsMV0sIFsxLDBdLCBbMCwyXSwgWzIsMF0sIFsxLDJdLCBbMiwxXV07XG52YXIgbmJFZGdlcyA9IDA7XG5jZWxscy5mb3JFYWNoKGZ1bmN0aW9uKGNlbGwpe1xuICBmb3IodmFyIGk9MDsgaTxwZXJtdXRhdGlvbnMubGVuZ3RoOyArK2kpe1xuICAgIHZhciBzID0gcGVybXV0YXRpb25zW2ldWzBdO1xuICAgIHZhciBkID0gcGVybXV0YXRpb25zW2ldWzFdO1xuICAgIHZhciBwcyA9IHBvaW50c1tjZWxsW3NdXTtcbiAgICB2YXIgcGQgPSBwb2ludHNbY2VsbFtkXV07XG4gICAgdmFyIGRpc3RhbmNlID0gc3FydCggcG93KHBzLnggLSBwZC54LCAyKSArIHBvdyhwcy55IC0gcGQueSwgMikgKTtcbiAgICB2YXIgZWRnZSA9IHtpZCA6IG5iRWRnZXMsXG4gICAgICAgICAgc291cmNlOiBjZWxsW3NdLCBcbiAgICAgICAgICBkZXN0aW5hdGlvbjogY2VsbFtkXSwgXG4gICAgICAgICAgZGlzdGFuY2UgOiBkaXN0YW5jZSxcbiAgICAgICAgICBkaXJlY3Rpb24gOiBNYXRoLmF0YW4oKHBkLnktcHMueSkvKHBkLngtcHMueCkpLFxuICAgICAgICAgIG9yaWVudGF0aW9uIDogc2lnbigocGQueC1wcy54KSksXG4gICAgICAgICAgcGhlcm9tb24gOiAxL2Rpc3RhbmNlXG4gICAgICAgICAgfTtcbiAgICB2YXIgbmV4dHM7XG4gICAgaWYobmV4dEVkZ2VzLmhhcyhwcy5pZCkpe1xuICAgICAgICBuZXh0cyA9IG5leHRFZGdlcy5nZXQocHMuaWQpO1xuICAgICAgICBuZXh0cy5wdXNoKGVkZ2UpO1xuICAgIH0gZWxzZSB7XG4gICAgICAgIG5leHRzID0gW2VkZ2VdO1xuICAgIH1cbiAgICBuZXh0RWRnZXMuc2V0KHBzLmlkLCBuZXh0cyk7XG4gICAgZWRnZXMucHVzaChlZGdlKTtcbiAgICBuYkVkZ2VzKys7XG4gIH1cbiAgXG59KVxuXG4vLyBpbml0aWFsaXplIGFudHNcbnZhciBwb3B1bGF0aW9uID0gbmV3IEFycmF5KG5iQW50cyk7XG52YXIgaSxqO1xuZm9yIChpID0gMDsgaSA8IG5iQW50czsgaSsrKSB7XG4gICAgLy8gdGFrZSBhIHJhbmRvbSBlZGdlXG4gICAgdmFyIGVkZ2UgPSBlZGdlc1tNYXRoLmZsb29yKGVkZ2VzLmxlbmd0aCpyYW5kb20oKSldO1xuICAgIHZhciB4ID0gcG9pbnRzW2VkZ2Uuc291cmNlXS54IFxuICAgIHZhciB5ID0gcG9pbnRzW2VkZ2Uuc291cmNlXS55XG4gICAgcG9wdWxhdGlvbltpXSA9IG5ldyBBbnQoeCwgeSwgZWRnZSk7XG59XG5cbnZhciBjYW52YXMsIGNvbnRleHQ7XG5cbnNoZWxsLm9uKFwiaW5pdFwiLCBmdW5jdGlvbigpIHtcbiAgICBjYW52YXMgPSBkb2N1bWVudC5jcmVhdGVFbGVtZW50KFwiY2FudmFzXCIpO1xuICAgIGNhbnZhcy53aWR0aCA9IHNoZWxsLndpZHRoO1xuICAgIGNhbnZhcy5oZWlnaHQgPSBzaGVsbC5oZWlnaHQ7XG4gICAgY29udGV4dCA9IGNhbnZhcy5nZXRDb250ZXh0KFwiMmRcIik7XG4gICAgc2hlbGwuZWxlbWVudC5hcHBlbmRDaGlsZChjYW52YXMpO1xufSlcblxuc2hlbGwub24oXCJyZXNpemVcIiwgZnVuY3Rpb24odywgaCkge1xuICAgIGNhbnZhcy53aWR0aCA9IHc7XG4gICAgY2FudmFzLmhlaWdodCA9IGg7XG59KVxuXG5zaGVsbC5vbihcInJlbmRlclwiLCBmdW5jdGlvbigpIHtcbiAgICB2YXIgdyA9IGNhbnZhcy53aWR0aDtcbiAgICB2YXIgaCA9IGNhbnZhcy5oZWlnaHQ7XG4gICAgdmFyIG1vdXNlID0gW3NoZWxsLm1vdXNlWC93LCBzaGVsbC5tb3VzZVkvaF07XG4gICAgY29udGV4dC5zZXRUcmFuc2Zvcm0odywgMCwgMCwgaCwgMCwgMCk7XG4gICAgY29udGV4dC5maWxsU3R5bGUgPSBcIiNmZmZcIjtcbiAgICBjb250ZXh0LmZpbGxSZWN0KDAsMCx3LGgpO1xuXG4gICAgLy8gZWRnZXNcbiAgICBjb250ZXh0LnN0cm9rZVN0eWxlID0gXCIjMDAwXCI7XG4gICAgZm9yKHZhciBpPTA7IGk8ZWRnZXMubGVuZ3RoOyArK2kpIHtcbiAgICAgICAgdmFyIGVkZ2UgPSBlZGdlc1tpXTtcbiAgICAgICAgaWYgKGVkZ2UucGhlcm9tb24gIT0gMCl7XG4gICAgICAgICAgICBjb250ZXh0LmxpbmVXaWR0aCA9IDAuMDAwMDEgKiBlZGdlLnBoZXJvbW9uO1xuICAgICAgICB9ZWxzZSB7XG4gICAgICAgICAgICBjb250ZXh0LmxpbmVXaWR0aCA9IDAuMDAwMDE7XG4gICAgICAgIH1cbiAgICAgICAgY29udGV4dC5iZWdpblBhdGgoKTtcbiAgICAgICAgY29udGV4dC5tb3ZlVG8ocG9pbnRzW2VkZ2Uuc291cmNlXS54LCBwb2ludHNbZWRnZS5zb3VyY2VdLnkpO1xuICAgICAgICBjb250ZXh0LmxpbmVUbyhwb2ludHNbZWRnZS5kZXN0aW5hdGlvbl0ueCwgcG9pbnRzW2VkZ2UuZGVzdGluYXRpb25dLnkpO1xuICAgICAgICBjb250ZXh0LnN0cm9rZSgpO1xuICAgIH1cblxuICAgIC8vIHZlcnRpY2VzXG4gICAgZm9yKHZhciBpPTA7IGk8cG9pbnRzLmxlbmd0aDsgKytpKSB7XG4gICAgICAgIGNvbnRleHQuYmVnaW5QYXRoKClcbiAgICAgICAgdmFyIHBvaW50ID0gcG9pbnRzW2ldO1xuICAgICAgICBpZiAoY2l0eVNldC5oYXMocG9pbnQuaWQpKSB7XG4gICAgICAgICAgICBjb250ZXh0LmZpbGxTdHlsZSA9IFwiIzAxMDFERlwiO1xuICAgICAgICAgICAgY29udGV4dC5hcmMocG9pbnQueCwgcG9pbnQueSwgMC4wMDYsIDAsIDIqTWF0aC5QSSk7XG4gICAgICAgIH1cbiAgICAgICAgZWxzZSB7XG4gICAgICAgICAgICBjb250ZXh0LmZpbGxTdHlsZSA9IFwiIzAwMFwiO1xuICAgICAgICAgICAgY29udGV4dC5hcmMocG9pbnRzW2ldLngsIHBvaW50c1tpXS55LCAwLjAwMywgMCwgMipNYXRoLlBJKTtcbiAgICAgICAgfVxuICAgICAgICBjb250ZXh0LmNsb3NlUGF0aCgpO1xuICAgICAgICBjb250ZXh0LmZpbGwoKTtcbiAgICB9XG5cbiAgICAvLyBtb3ZlIGFudHNcbiAgICBmb3IgKGkgPSAwOyBpIDwgbmJBbnRzOyBpKyspIHtcbiAgICAgICAgcG9wdWxhdGlvbltpXS50cmFuc2l0KCk7XG4gICAgfVxuXG4gICAgLy8gcGhlcm9tb24gZXZhcG9yYXRpb25cbiAgICBmb3IgKGkgPSAwOyBpIDwgZWRnZXMubGVuZ3RoOyBpKyspIHtcbiAgICAgICAgaWYoZWRnZXNbaV0ucGhlcm9tb24gPiAwKXtcbiAgICAgICAgICAgIGVkZ2VzW2ldLnBoZXJvbW9uIC09IDAuMDAxO1xuICAgICAgICB9XG4gICAgfVxuXG5cbiAgICBmb3IodmFyIGk9MDsgaTxwb3B1bGF0aW9uLmxlbmd0aDsgKytpKSB7XG4gICAgICAgIGNvbnRleHQuYmVnaW5QYXRoKClcbiAgICAgICAgdmFyIHggPSBwb3B1bGF0aW9uW2ldLnBvc1ggLy8rIDAuMDEqcmFuZG9tKCk7XG4gICAgICAgIHZhciB5ID0gcG9wdWxhdGlvbltpXS5wb3NZIC8vKyAwLjAxKnJhbmRvbSgpO1xuICAgICAgICBpZiAocG9wdWxhdGlvbltpXS5zdGF0ZSA9PT0gXCJwaGVyb21vbmluZ1wiKXtjb250ZXh0LmZpbGxTdHlsZSA9IFwiI0ZGMDAwMFwifVxuICAgICAgICBlbHNlIHtjb250ZXh0LmZpbGxTdHlsZSA9IFwiIzYxMEIwQlwifVxuICAgICAgICBjb250ZXh0LmFyYyh4LCB5LCAwLjAwMywgMCwgMipNYXRoLlBJKVxuICAgICAgICBjb250ZXh0LmNsb3NlUGF0aCgpXG4gICAgICAgIGNvbnRleHQuZmlsbCgpXG4gICAgfVxuICBcbn0pXG5cbmZ1bmN0aW9uIEFudCh4LCB5LCBlZGdlKSB7ICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICBcbiAgICB0aGlzLnBvc1ggPSB4OyAgICAgICAgICAgICAgICBcbiAgICB0aGlzLnBvc1kgPSB5O1xuICAgIHRoaXMuZWRnZSA9IGVkZ2U7XG4gICAgdGhpcy5zdGVwID0gMDtcbiAgICB0aGlzLnN0YXRlID0gXCJmb3JhZ2VcIjtcbiAgICB0aGlzLnRyYW5zaXQgPSBzdGF0ZW1hY2hpbmU7IFxuICAgIHRoaXMubW92ZSA9IG1vdmU7XG4gICAgdGhpcy5lZGdlcyA9IFtdO1xuICAgIHRoaXMubGFzdENpdHkgPSB1bmRlZmluZWQ7XG59XG4vLyBmb3JhZ2U6IHRoZSBhbnQgd2FuZGVycyBhcm91bmQgd2l0aG91dCBhbnkgcGhlcm9tb24gZGVwb3NpdGlvblxuLy8gb25jZSBpdCBmaW5kcyBhIGNpdHksIGl0IHN0YXJ0cyByZW1lbWJlcmluZyB0aGUgbm9kZXMgaXQgZ29lcyB0aHJvdWdoXG4vLyB3aGVuIGl0IGZpbmRzIGFub3RoZXIgY2l0eSwgaXQgY29tcHV0ZXMgdGhlIHBhdGggbGVuZ3RoIGFuZCBhZGRzIHBoZXJvbW9ucyBvbmUgZWFjaCBlZGdlc1xuLy8gcHJvcG9ydGlvbm5hbHkgdG8gdGhlIHNob3J0ZXN0bmVzcyBvZiB0aGUgcGF0aFxuLy8gaXQgcmVzZXRzIHRoZSBsaXN0IG9mIG5vZGVzIGFuZCBjb250aW51ZXNcbi8vIHdoaWxlIGZvcmFnaW5nIHRoZSBhbnQgY2hvc2VzIHRoZSBwYXRoIHdpdGggYSBwaGVyb21vbiBwcmVmZXJlbmNlXG5cbmZ1bmN0aW9uIG1vdmUoKSB7XG4gICAgdmFyIGVkZ2VDaGFuZ2VkO1xuICAgIHZhciBjaXR5UmVhY2hlZCA9IGZhbHNlO1xuICAgIC8vIG9uIGVkZ2VcbiAgICBpZiAodGhpcy5zdGVwIDwgdGhpcy5lZGdlLmRpc3RhbmNlKXtcbiAgICAgICAgdGhpcy5wb3NYICs9IDAuMDA1Kk1hdGguY29zKHRoaXMuZWRnZS5kaXJlY3Rpb24pKnRoaXMuZWRnZS5vcmllbnRhdGlvbjtcbiAgICAgICAgdGhpcy5wb3NZICs9IDAuMDA1Kk1hdGguc2luKHRoaXMuZWRnZS5kaXJlY3Rpb24pKnRoaXMuZWRnZS5vcmllbnRhdGlvbjtcbiAgICAgICAgdGhpcy5zdGVwICs9IDAuMDA1O1xuICAgICAgICBlZGdlQ2hhbmdlZCA9IGZhbHNlO1xuICAgIC8vIG9uIHZlcnRleFxuICAgIH0gZWxzZSB7XG4gICAgICAgIHRoaXMuc3RlcCA9IDA7XG4gICAgICAgIHRoaXMucG9zWCA9IHBvaW50c1t0aGlzLmVkZ2UuZGVzdGluYXRpb25dLng7XG4gICAgICAgIHRoaXMucG9zWSA9IHBvaW50c1t0aGlzLmVkZ2UuZGVzdGluYXRpb25dLnk7XG4gICAgICAgIHZhciBwb3NzaWJsZUVkZ2VzID0gbmV4dEVkZ2VzLmdldCh0aGlzLmVkZ2UuZGVzdGluYXRpb24pO1xuICAgICAgICAvLyBmbGlwIGEgY29pbiBhbmQgZWl0aGVyIHRha2UgdGhlIHNtZWxsaWVzdCBwYXRoIG9mIGEgcmFuZG9tIG9uZVxuICAgICAgICBpZiAocmFuZG9tKCkgPiAwLjUpe1xuICAgICAgICAgICAgdmFyIHNtZWxscyA9IHBvc3NpYmxlRWRnZXMubWFwKGZ1bmN0aW9uKGUpe3JldHVybiBlLnBoZXJvbW9ufSk7XG4gICAgICAgICAgICB2YXIgaW5kZXggPSBzbWVsbHMuaW5kZXhPZihNYXRoLm1heC5hcHBseShNYXRoLCBzbWVsbHMpKTtcbiAgICAgICAgICAgIHRoaXMuZWRnZSA9IHBvc3NpYmxlRWRnZXNbaW5kZXhdO1xuICAgICAgICB9IGVsc2Uge1xuICAgICAgICAgICAgdGhpcy5lZGdlID0gcG9zc2libGVFZGdlc1tmbG9vcihyYW5kb20oKSpwb3NzaWJsZUVkZ2VzLmxlbmd0aCldO1xuICAgICAgICB9XG4gICAgICAgIGNpdHlSZWFjaGVkID0gY2l0eVNldC5oYXModGhpcy5lZGdlLnNvdXJjZSk7XG4gICAgICAgIGVkZ2VDaGFuZ2VkID0gdHJ1ZTtcbiAgICB9XG4gICAgcmV0dXJuIHtjaXR5UmVhY2hlZDogY2l0eVJlYWNoZWQsIGVkZ2VDaGFuZ2VkOiBlZGdlQ2hhbmdlZH07XG59XG5cblxuZnVuY3Rpb24gc3RhdGVtYWNoaW5lKCkge1xuICAgIHN3aXRjaCAodGhpcy5zdGF0ZSkge1xuICAgICAgICBjYXNlIFwiZm9yYWdlXCI6XG4gICAgICAgICAgICB2YXIgcmVzID0gdGhpcy5tb3ZlKCk7XG4gICAgICAgICAgICBpZiAocmVzLmNpdHlSZWFjaGVkKSB7XG4gICAgICAgICAgICAgICAgdGhpcy5zdGF0ZSA9IFwicGhlcm9tb25pbmdcIjtcbiAgICAgICAgICAgICAgICB0aGlzLmxhc3RDaXR5ID0gdGhpcy5lZGdlLnNvdXJjZTtcbiAgICAgICAgICAgIH07XG4gICAgICAgICAgICBicmVhaztcbiAgICAgICAgY2FzZSBcInBoZXJvbW9uaW5nXCI6XG4gICAgICAgICAgICB2YXIgcmVzID0gdGhpcy5tb3ZlKCk7XG4gICAgICAgICAgICBpZiAocmVzLmVkZ2VDaGFuZ2VkKSB7XG4gICAgICAgICAgICAgICAgdGhpcy5lZGdlcy5wdXNoKHRoaXMuZWRnZSk7XG4gICAgICAgICAgICAgICAgLy8gZm91bmQgYSBjaXR5XG4gICAgICAgICAgICAgICAgaWYgKHJlcy5jaXR5UmVhY2hlZCAmJiAodGhpcy5lZGdlLnNvdXJjZSAhPSB0aGlzLmxhc3RDaXR5KSApe1xuICAgICAgICAgICAgICAgICAgICAvLyBjb21wdXRlIHRoZSBsZW5ndGggb2YgdGhlIHBhdGhcbiAgICAgICAgICAgICAgICAgICAgdmFyIHBhdGhMZW5ndGggPSB0aGlzLmVkZ2VzLm1hcChmdW5jdGlvbihlKXtyZXR1cm4gZS5kaXN0YW5jZX0pLnJlZHVjZShmdW5jdGlvbihhLGIpe3JldHVybiBhICsgYn0pO1xuICAgICAgICAgICAgICAgICAgICB2YXIgZGVsdGFQaGVyb21vbmUgPSAxL3BhdGhMZW5ndGg7XG4gICAgICAgICAgICAgICAgICAgIHRoaXMuZWRnZXMuZm9yRWFjaChmdW5jdGlvbihlKXtlLnBoZXJvbW9uICs9IGRlbHRhUGhlcm9tb25lfSk7XG4gICAgICAgICAgICAgICAgICAgIC8vIGNvbnNvbGUubG9nKGRlbHRhUGhlcm9tb25lLCB0aGlzLmVkZ2VzKTtcbiAgICAgICAgICAgICAgICAgICAgdGhpcy5lZGdlcyA9IFt0aGlzLmVkZ2VdO1xuICAgICAgICAgICAgICAgICAgICB0aGlzLmxhc3RDaXR5ID0gdGhpcy5lZGdlLnNvdXJjZTtcbiAgICAgICAgICAgICAgICB9XG4gICAgICAgICAgICB9XG4gICAgICAgICAgYnJlYWs7XG5cbiAgICB9XG59IiwiLy8gQ29weXJpZ2h0IEpveWVudCwgSW5jLiBhbmQgb3RoZXIgTm9kZSBjb250cmlidXRvcnMuXG4vL1xuLy8gUGVybWlzc2lvbiBpcyBoZXJlYnkgZ3JhbnRlZCwgZnJlZSBvZiBjaGFyZ2UsIHRvIGFueSBwZXJzb24gb2J0YWluaW5nIGFcbi8vIGNvcHkgb2YgdGhpcyBzb2Z0d2FyZSBhbmQgYXNzb2NpYXRlZCBkb2N1bWVudGF0aW9uIGZpbGVzICh0aGVcbi8vIFwiU29mdHdhcmVcIiksIHRvIGRlYWwgaW4gdGhlIFNvZnR3YXJlIHdpdGhvdXQgcmVzdHJpY3Rpb24sIGluY2x1ZGluZ1xuLy8gd2l0aG91dCBsaW1pdGF0aW9uIHRoZSByaWdodHMgdG8gdXNlLCBjb3B5LCBtb2RpZnksIG1lcmdlLCBwdWJsaXNoLFxuLy8gZGlzdHJpYnV0ZSwgc3VibGljZW5zZSwgYW5kL29yIHNlbGwgY29waWVzIG9mIHRoZSBTb2Z0d2FyZSwgYW5kIHRvIHBlcm1pdFxuLy8gcGVyc29ucyB0byB3aG9tIHRoZSBTb2Z0d2FyZSBpcyBmdXJuaXNoZWQgdG8gZG8gc28sIHN1YmplY3QgdG8gdGhlXG4vLyBmb2xsb3dpbmcgY29uZGl0aW9uczpcbi8vXG4vLyBUaGUgYWJvdmUgY29weXJpZ2h0IG5vdGljZSBhbmQgdGhpcyBwZXJtaXNzaW9uIG5vdGljZSBzaGFsbCBiZSBpbmNsdWRlZFxuLy8gaW4gYWxsIGNvcGllcyBvciBzdWJzdGFudGlhbCBwb3J0aW9ucyBvZiB0aGUgU29mdHdhcmUuXG4vL1xuLy8gVEhFIFNPRlRXQVJFIElTIFBST1ZJREVEIFwiQVMgSVNcIiwgV0lUSE9VVCBXQVJSQU5UWSBPRiBBTlkgS0lORCwgRVhQUkVTU1xuLy8gT1IgSU1QTElFRCwgSU5DTFVESU5HIEJVVCBOT1QgTElNSVRFRCBUTyBUSEUgV0FSUkFOVElFUyBPRlxuLy8gTUVSQ0hBTlRBQklMSVRZLCBGSVRORVNTIEZPUiBBIFBBUlRJQ1VMQVIgUFVSUE9TRSBBTkQgTk9OSU5GUklOR0VNRU5ULiBJTlxuLy8gTk8gRVZFTlQgU0hBTEwgVEhFIEFVVEhPUlMgT1IgQ09QWVJJR0hUIEhPTERFUlMgQkUgTElBQkxFIEZPUiBBTlkgQ0xBSU0sXG4vLyBEQU1BR0VTIE9SIE9USEVSIExJQUJJTElUWSwgV0hFVEhFUiBJTiBBTiBBQ1RJT04gT0YgQ09OVFJBQ1QsIFRPUlQgT1Jcbi8vIE9USEVSV0lTRSwgQVJJU0lORyBGUk9NLCBPVVQgT0YgT1IgSU4gQ09OTkVDVElPTiBXSVRIIFRIRSBTT0ZUV0FSRSBPUiBUSEVcbi8vIFVTRSBPUiBPVEhFUiBERUFMSU5HUyBJTiBUSEUgU09GVFdBUkUuXG5cbmZ1bmN0aW9uIEV2ZW50RW1pdHRlcigpIHtcbiAgdGhpcy5fZXZlbnRzID0gdGhpcy5fZXZlbnRzIHx8IHt9O1xuICB0aGlzLl9tYXhMaXN0ZW5lcnMgPSB0aGlzLl9tYXhMaXN0ZW5lcnMgfHwgdW5kZWZpbmVkO1xufVxubW9kdWxlLmV4cG9ydHMgPSBFdmVudEVtaXR0ZXI7XG5cbi8vIEJhY2t3YXJkcy1jb21wYXQgd2l0aCBub2RlIDAuMTAueFxuRXZlbnRFbWl0dGVyLkV2ZW50RW1pdHRlciA9IEV2ZW50RW1pdHRlcjtcblxuRXZlbnRFbWl0dGVyLnByb3RvdHlwZS5fZXZlbnRzID0gdW5kZWZpbmVkO1xuRXZlbnRFbWl0dGVyLnByb3RvdHlwZS5fbWF4TGlzdGVuZXJzID0gdW5kZWZpbmVkO1xuXG4vLyBCeSBkZWZhdWx0IEV2ZW50RW1pdHRlcnMgd2lsbCBwcmludCBhIHdhcm5pbmcgaWYgbW9yZSB0aGFuIDEwIGxpc3RlbmVycyBhcmVcbi8vIGFkZGVkIHRvIGl0LiBUaGlzIGlzIGEgdXNlZnVsIGRlZmF1bHQgd2hpY2ggaGVscHMgZmluZGluZyBtZW1vcnkgbGVha3MuXG5FdmVudEVtaXR0ZXIuZGVmYXVsdE1heExpc3RlbmVycyA9IDEwO1xuXG4vLyBPYnZpb3VzbHkgbm90IGFsbCBFbWl0dGVycyBzaG91bGQgYmUgbGltaXRlZCB0byAxMC4gVGhpcyBmdW5jdGlvbiBhbGxvd3Ncbi8vIHRoYXQgdG8gYmUgaW5jcmVhc2VkLiBTZXQgdG8gemVybyBmb3IgdW5saW1pdGVkLlxuRXZlbnRFbWl0dGVyLnByb3RvdHlwZS5zZXRNYXhMaXN0ZW5lcnMgPSBmdW5jdGlvbihuKSB7XG4gIGlmICghaXNOdW1iZXIobikgfHwgbiA8IDAgfHwgaXNOYU4obikpXG4gICAgdGhyb3cgVHlwZUVycm9yKCduIG11c3QgYmUgYSBwb3NpdGl2ZSBudW1iZXInKTtcbiAgdGhpcy5fbWF4TGlzdGVuZXJzID0gbjtcbiAgcmV0dXJuIHRoaXM7XG59O1xuXG5FdmVudEVtaXR0ZXIucHJvdG90eXBlLmVtaXQgPSBmdW5jdGlvbih0eXBlKSB7XG4gIHZhciBlciwgaGFuZGxlciwgbGVuLCBhcmdzLCBpLCBsaXN0ZW5lcnM7XG5cbiAgaWYgKCF0aGlzLl9ldmVudHMpXG4gICAgdGhpcy5fZXZlbnRzID0ge307XG5cbiAgLy8gSWYgdGhlcmUgaXMgbm8gJ2Vycm9yJyBldmVudCBsaXN0ZW5lciB0aGVuIHRocm93LlxuICBpZiAodHlwZSA9PT0gJ2Vycm9yJykge1xuICAgIGlmICghdGhpcy5fZXZlbnRzLmVycm9yIHx8XG4gICAgICAgIChpc09iamVjdCh0aGlzLl9ldmVudHMuZXJyb3IpICYmICF0aGlzLl9ldmVudHMuZXJyb3IubGVuZ3RoKSkge1xuICAgICAgZXIgPSBhcmd1bWVudHNbMV07XG4gICAgICBpZiAoZXIgaW5zdGFuY2VvZiBFcnJvcikge1xuICAgICAgICB0aHJvdyBlcjsgLy8gVW5oYW5kbGVkICdlcnJvcicgZXZlbnRcbiAgICAgIH0gZWxzZSB7XG4gICAgICAgIHRocm93IFR5cGVFcnJvcignVW5jYXVnaHQsIHVuc3BlY2lmaWVkIFwiZXJyb3JcIiBldmVudC4nKTtcbiAgICAgIH1cbiAgICAgIHJldHVybiBmYWxzZTtcbiAgICB9XG4gIH1cblxuICBoYW5kbGVyID0gdGhpcy5fZXZlbnRzW3R5cGVdO1xuXG4gIGlmIChpc1VuZGVmaW5lZChoYW5kbGVyKSlcbiAgICByZXR1cm4gZmFsc2U7XG5cbiAgaWYgKGlzRnVuY3Rpb24oaGFuZGxlcikpIHtcbiAgICBzd2l0Y2ggKGFyZ3VtZW50cy5sZW5ndGgpIHtcbiAgICAgIC8vIGZhc3QgY2FzZXNcbiAgICAgIGNhc2UgMTpcbiAgICAgICAgaGFuZGxlci5jYWxsKHRoaXMpO1xuICAgICAgICBicmVhaztcbiAgICAgIGNhc2UgMjpcbiAgICAgICAgaGFuZGxlci5jYWxsKHRoaXMsIGFyZ3VtZW50c1sxXSk7XG4gICAgICAgIGJyZWFrO1xuICAgICAgY2FzZSAzOlxuICAgICAgICBoYW5kbGVyLmNhbGwodGhpcywgYXJndW1lbnRzWzFdLCBhcmd1bWVudHNbMl0pO1xuICAgICAgICBicmVhaztcbiAgICAgIC8vIHNsb3dlclxuICAgICAgZGVmYXVsdDpcbiAgICAgICAgbGVuID0gYXJndW1lbnRzLmxlbmd0aDtcbiAgICAgICAgYXJncyA9IG5ldyBBcnJheShsZW4gLSAxKTtcbiAgICAgICAgZm9yIChpID0gMTsgaSA8IGxlbjsgaSsrKVxuICAgICAgICAgIGFyZ3NbaSAtIDFdID0gYXJndW1lbnRzW2ldO1xuICAgICAgICBoYW5kbGVyLmFwcGx5KHRoaXMsIGFyZ3MpO1xuICAgIH1cbiAgfSBlbHNlIGlmIChpc09iamVjdChoYW5kbGVyKSkge1xuICAgIGxlbiA9IGFyZ3VtZW50cy5sZW5ndGg7XG4gICAgYXJncyA9IG5ldyBBcnJheShsZW4gLSAxKTtcbiAgICBmb3IgKGkgPSAxOyBpIDwgbGVuOyBpKyspXG4gICAgICBhcmdzW2kgLSAxXSA9IGFyZ3VtZW50c1tpXTtcblxuICAgIGxpc3RlbmVycyA9IGhhbmRsZXIuc2xpY2UoKTtcbiAgICBsZW4gPSBsaXN0ZW5lcnMubGVuZ3RoO1xuICAgIGZvciAoaSA9IDA7IGkgPCBsZW47IGkrKylcbiAgICAgIGxpc3RlbmVyc1tpXS5hcHBseSh0aGlzLCBhcmdzKTtcbiAgfVxuXG4gIHJldHVybiB0cnVlO1xufTtcblxuRXZlbnRFbWl0dGVyLnByb3RvdHlwZS5hZGRMaXN0ZW5lciA9IGZ1bmN0aW9uKHR5cGUsIGxpc3RlbmVyKSB7XG4gIHZhciBtO1xuXG4gIGlmICghaXNGdW5jdGlvbihsaXN0ZW5lcikpXG4gICAgdGhyb3cgVHlwZUVycm9yKCdsaXN0ZW5lciBtdXN0IGJlIGEgZnVuY3Rpb24nKTtcblxuICBpZiAoIXRoaXMuX2V2ZW50cylcbiAgICB0aGlzLl9ldmVudHMgPSB7fTtcblxuICAvLyBUbyBhdm9pZCByZWN1cnNpb24gaW4gdGhlIGNhc2UgdGhhdCB0eXBlID09PSBcIm5ld0xpc3RlbmVyXCIhIEJlZm9yZVxuICAvLyBhZGRpbmcgaXQgdG8gdGhlIGxpc3RlbmVycywgZmlyc3QgZW1pdCBcIm5ld0xpc3RlbmVyXCIuXG4gIGlmICh0aGlzLl9ldmVudHMubmV3TGlzdGVuZXIpXG4gICAgdGhpcy5lbWl0KCduZXdMaXN0ZW5lcicsIHR5cGUsXG4gICAgICAgICAgICAgIGlzRnVuY3Rpb24obGlzdGVuZXIubGlzdGVuZXIpID9cbiAgICAgICAgICAgICAgbGlzdGVuZXIubGlzdGVuZXIgOiBsaXN0ZW5lcik7XG5cbiAgaWYgKCF0aGlzLl9ldmVudHNbdHlwZV0pXG4gICAgLy8gT3B0aW1pemUgdGhlIGNhc2Ugb2Ygb25lIGxpc3RlbmVyLiBEb24ndCBuZWVkIHRoZSBleHRyYSBhcnJheSBvYmplY3QuXG4gICAgdGhpcy5fZXZlbnRzW3R5cGVdID0gbGlzdGVuZXI7XG4gIGVsc2UgaWYgKGlzT2JqZWN0KHRoaXMuX2V2ZW50c1t0eXBlXSkpXG4gICAgLy8gSWYgd2UndmUgYWxyZWFkeSBnb3QgYW4gYXJyYXksIGp1c3QgYXBwZW5kLlxuICAgIHRoaXMuX2V2ZW50c1t0eXBlXS5wdXNoKGxpc3RlbmVyKTtcbiAgZWxzZVxuICAgIC8vIEFkZGluZyB0aGUgc2Vjb25kIGVsZW1lbnQsIG5lZWQgdG8gY2hhbmdlIHRvIGFycmF5LlxuICAgIHRoaXMuX2V2ZW50c1t0eXBlXSA9IFt0aGlzLl9ldmVudHNbdHlwZV0sIGxpc3RlbmVyXTtcblxuICAvLyBDaGVjayBmb3IgbGlzdGVuZXIgbGVha1xuICBpZiAoaXNPYmplY3QodGhpcy5fZXZlbnRzW3R5cGVdKSAmJiAhdGhpcy5fZXZlbnRzW3R5cGVdLndhcm5lZCkge1xuICAgIHZhciBtO1xuICAgIGlmICghaXNVbmRlZmluZWQodGhpcy5fbWF4TGlzdGVuZXJzKSkge1xuICAgICAgbSA9IHRoaXMuX21heExpc3RlbmVycztcbiAgICB9IGVsc2Uge1xuICAgICAgbSA9IEV2ZW50RW1pdHRlci5kZWZhdWx0TWF4TGlzdGVuZXJzO1xuICAgIH1cblxuICAgIGlmIChtICYmIG0gPiAwICYmIHRoaXMuX2V2ZW50c1t0eXBlXS5sZW5ndGggPiBtKSB7XG4gICAgICB0aGlzLl9ldmVudHNbdHlwZV0ud2FybmVkID0gdHJ1ZTtcbiAgICAgIGNvbnNvbGUuZXJyb3IoJyhub2RlKSB3YXJuaW5nOiBwb3NzaWJsZSBFdmVudEVtaXR0ZXIgbWVtb3J5ICcgK1xuICAgICAgICAgICAgICAgICAgICAnbGVhayBkZXRlY3RlZC4gJWQgbGlzdGVuZXJzIGFkZGVkLiAnICtcbiAgICAgICAgICAgICAgICAgICAgJ1VzZSBlbWl0dGVyLnNldE1heExpc3RlbmVycygpIHRvIGluY3JlYXNlIGxpbWl0LicsXG4gICAgICAgICAgICAgICAgICAgIHRoaXMuX2V2ZW50c1t0eXBlXS5sZW5ndGgpO1xuICAgICAgaWYgKHR5cGVvZiBjb25zb2xlLnRyYWNlID09PSAnZnVuY3Rpb24nKSB7XG4gICAgICAgIC8vIG5vdCBzdXBwb3J0ZWQgaW4gSUUgMTBcbiAgICAgICAgY29uc29sZS50cmFjZSgpO1xuICAgICAgfVxuICAgIH1cbiAgfVxuXG4gIHJldHVybiB0aGlzO1xufTtcblxuRXZlbnRFbWl0dGVyLnByb3RvdHlwZS5vbiA9IEV2ZW50RW1pdHRlci5wcm90b3R5cGUuYWRkTGlzdGVuZXI7XG5cbkV2ZW50RW1pdHRlci5wcm90b3R5cGUub25jZSA9IGZ1bmN0aW9uKHR5cGUsIGxpc3RlbmVyKSB7XG4gIGlmICghaXNGdW5jdGlvbihsaXN0ZW5lcikpXG4gICAgdGhyb3cgVHlwZUVycm9yKCdsaXN0ZW5lciBtdXN0IGJlIGEgZnVuY3Rpb24nKTtcblxuICB2YXIgZmlyZWQgPSBmYWxzZTtcblxuICBmdW5jdGlvbiBnKCkge1xuICAgIHRoaXMucmVtb3ZlTGlzdGVuZXIodHlwZSwgZyk7XG5cbiAgICBpZiAoIWZpcmVkKSB7XG4gICAgICBmaXJlZCA9IHRydWU7XG4gICAgICBsaXN0ZW5lci5hcHBseSh0aGlzLCBhcmd1bWVudHMpO1xuICAgIH1cbiAgfVxuXG4gIGcubGlzdGVuZXIgPSBsaXN0ZW5lcjtcbiAgdGhpcy5vbih0eXBlLCBnKTtcblxuICByZXR1cm4gdGhpcztcbn07XG5cbi8vIGVtaXRzIGEgJ3JlbW92ZUxpc3RlbmVyJyBldmVudCBpZmYgdGhlIGxpc3RlbmVyIHdhcyByZW1vdmVkXG5FdmVudEVtaXR0ZXIucHJvdG90eXBlLnJlbW92ZUxpc3RlbmVyID0gZnVuY3Rpb24odHlwZSwgbGlzdGVuZXIpIHtcbiAgdmFyIGxpc3QsIHBvc2l0aW9uLCBsZW5ndGgsIGk7XG5cbiAgaWYgKCFpc0Z1bmN0aW9uKGxpc3RlbmVyKSlcbiAgICB0aHJvdyBUeXBlRXJyb3IoJ2xpc3RlbmVyIG11c3QgYmUgYSBmdW5jdGlvbicpO1xuXG4gIGlmICghdGhpcy5fZXZlbnRzIHx8ICF0aGlzLl9ldmVudHNbdHlwZV0pXG4gICAgcmV0dXJuIHRoaXM7XG5cbiAgbGlzdCA9IHRoaXMuX2V2ZW50c1t0eXBlXTtcbiAgbGVuZ3RoID0gbGlzdC5sZW5ndGg7XG4gIHBvc2l0aW9uID0gLTE7XG5cbiAgaWYgKGxpc3QgPT09IGxpc3RlbmVyIHx8XG4gICAgICAoaXNGdW5jdGlvbihsaXN0Lmxpc3RlbmVyKSAmJiBsaXN0Lmxpc3RlbmVyID09PSBsaXN0ZW5lcikpIHtcbiAgICBkZWxldGUgdGhpcy5fZXZlbnRzW3R5cGVdO1xuICAgIGlmICh0aGlzLl9ldmVudHMucmVtb3ZlTGlzdGVuZXIpXG4gICAgICB0aGlzLmVtaXQoJ3JlbW92ZUxpc3RlbmVyJywgdHlwZSwgbGlzdGVuZXIpO1xuXG4gIH0gZWxzZSBpZiAoaXNPYmplY3QobGlzdCkpIHtcbiAgICBmb3IgKGkgPSBsZW5ndGg7IGktLSA+IDA7KSB7XG4gICAgICBpZiAobGlzdFtpXSA9PT0gbGlzdGVuZXIgfHxcbiAgICAgICAgICAobGlzdFtpXS5saXN0ZW5lciAmJiBsaXN0W2ldLmxpc3RlbmVyID09PSBsaXN0ZW5lcikpIHtcbiAgICAgICAgcG9zaXRpb24gPSBpO1xuICAgICAgICBicmVhaztcbiAgICAgIH1cbiAgICB9XG5cbiAgICBpZiAocG9zaXRpb24gPCAwKVxuICAgICAgcmV0dXJuIHRoaXM7XG5cbiAgICBpZiAobGlzdC5sZW5ndGggPT09IDEpIHtcbiAgICAgIGxpc3QubGVuZ3RoID0gMDtcbiAgICAgIGRlbGV0ZSB0aGlzLl9ldmVudHNbdHlwZV07XG4gICAgfSBlbHNlIHtcbiAgICAgIGxpc3Quc3BsaWNlKHBvc2l0aW9uLCAxKTtcbiAgICB9XG5cbiAgICBpZiAodGhpcy5fZXZlbnRzLnJlbW92ZUxpc3RlbmVyKVxuICAgICAgdGhpcy5lbWl0KCdyZW1vdmVMaXN0ZW5lcicsIHR5cGUsIGxpc3RlbmVyKTtcbiAgfVxuXG4gIHJldHVybiB0aGlzO1xufTtcblxuRXZlbnRFbWl0dGVyLnByb3RvdHlwZS5yZW1vdmVBbGxMaXN0ZW5lcnMgPSBmdW5jdGlvbih0eXBlKSB7XG4gIHZhciBrZXksIGxpc3RlbmVycztcblxuICBpZiAoIXRoaXMuX2V2ZW50cylcbiAgICByZXR1cm4gdGhpcztcblxuICAvLyBub3QgbGlzdGVuaW5nIGZvciByZW1vdmVMaXN0ZW5lciwgbm8gbmVlZCB0byBlbWl0XG4gIGlmICghdGhpcy5fZXZlbnRzLnJlbW92ZUxpc3RlbmVyKSB7XG4gICAgaWYgKGFyZ3VtZW50cy5sZW5ndGggPT09IDApXG4gICAgICB0aGlzLl9ldmVudHMgPSB7fTtcbiAgICBlbHNlIGlmICh0aGlzLl9ldmVudHNbdHlwZV0pXG4gICAgICBkZWxldGUgdGhpcy5fZXZlbnRzW3R5cGVdO1xuICAgIHJldHVybiB0aGlzO1xuICB9XG5cbiAgLy8gZW1pdCByZW1vdmVMaXN0ZW5lciBmb3IgYWxsIGxpc3RlbmVycyBvbiBhbGwgZXZlbnRzXG4gIGlmIChhcmd1bWVudHMubGVuZ3RoID09PSAwKSB7XG4gICAgZm9yIChrZXkgaW4gdGhpcy5fZXZlbnRzKSB7XG4gICAgICBpZiAoa2V5ID09PSAncmVtb3ZlTGlzdGVuZXInKSBjb250aW51ZTtcbiAgICAgIHRoaXMucmVtb3ZlQWxsTGlzdGVuZXJzKGtleSk7XG4gICAgfVxuICAgIHRoaXMucmVtb3ZlQWxsTGlzdGVuZXJzKCdyZW1vdmVMaXN0ZW5lcicpO1xuICAgIHRoaXMuX2V2ZW50cyA9IHt9O1xuICAgIHJldHVybiB0aGlzO1xuICB9XG5cbiAgbGlzdGVuZXJzID0gdGhpcy5fZXZlbnRzW3R5cGVdO1xuXG4gIGlmIChpc0Z1bmN0aW9uKGxpc3RlbmVycykpIHtcbiAgICB0aGlzLnJlbW92ZUxpc3RlbmVyKHR5cGUsIGxpc3RlbmVycyk7XG4gIH0gZWxzZSB7XG4gICAgLy8gTElGTyBvcmRlclxuICAgIHdoaWxlIChsaXN0ZW5lcnMubGVuZ3RoKVxuICAgICAgdGhpcy5yZW1vdmVMaXN0ZW5lcih0eXBlLCBsaXN0ZW5lcnNbbGlzdGVuZXJzLmxlbmd0aCAtIDFdKTtcbiAgfVxuICBkZWxldGUgdGhpcy5fZXZlbnRzW3R5cGVdO1xuXG4gIHJldHVybiB0aGlzO1xufTtcblxuRXZlbnRFbWl0dGVyLnByb3RvdHlwZS5saXN0ZW5lcnMgPSBmdW5jdGlvbih0eXBlKSB7XG4gIHZhciByZXQ7XG4gIGlmICghdGhpcy5fZXZlbnRzIHx8ICF0aGlzLl9ldmVudHNbdHlwZV0pXG4gICAgcmV0ID0gW107XG4gIGVsc2UgaWYgKGlzRnVuY3Rpb24odGhpcy5fZXZlbnRzW3R5cGVdKSlcbiAgICByZXQgPSBbdGhpcy5fZXZlbnRzW3R5cGVdXTtcbiAgZWxzZVxuICAgIHJldCA9IHRoaXMuX2V2ZW50c1t0eXBlXS5zbGljZSgpO1xuICByZXR1cm4gcmV0O1xufTtcblxuRXZlbnRFbWl0dGVyLmxpc3RlbmVyQ291bnQgPSBmdW5jdGlvbihlbWl0dGVyLCB0eXBlKSB7XG4gIHZhciByZXQ7XG4gIGlmICghZW1pdHRlci5fZXZlbnRzIHx8ICFlbWl0dGVyLl9ldmVudHNbdHlwZV0pXG4gICAgcmV0ID0gMDtcbiAgZWxzZSBpZiAoaXNGdW5jdGlvbihlbWl0dGVyLl9ldmVudHNbdHlwZV0pKVxuICAgIHJldCA9IDE7XG4gIGVsc2VcbiAgICByZXQgPSBlbWl0dGVyLl9ldmVudHNbdHlwZV0ubGVuZ3RoO1xuICByZXR1cm4gcmV0O1xufTtcblxuZnVuY3Rpb24gaXNGdW5jdGlvbihhcmcpIHtcbiAgcmV0dXJuIHR5cGVvZiBhcmcgPT09ICdmdW5jdGlvbic7XG59XG5cbmZ1bmN0aW9uIGlzTnVtYmVyKGFyZykge1xuICByZXR1cm4gdHlwZW9mIGFyZyA9PT0gJ251bWJlcic7XG59XG5cbmZ1bmN0aW9uIGlzT2JqZWN0KGFyZykge1xuICByZXR1cm4gdHlwZW9mIGFyZyA9PT0gJ29iamVjdCcgJiYgYXJnICE9PSBudWxsO1xufVxuXG5mdW5jdGlvbiBpc1VuZGVmaW5lZChhcmcpIHtcbiAgcmV0dXJuIGFyZyA9PT0gdm9pZCAwO1xufVxuIiwiaWYgKHR5cGVvZiBPYmplY3QuY3JlYXRlID09PSAnZnVuY3Rpb24nKSB7XG4gIC8vIGltcGxlbWVudGF0aW9uIGZyb20gc3RhbmRhcmQgbm9kZS5qcyAndXRpbCcgbW9kdWxlXG4gIG1vZHVsZS5leHBvcnRzID0gZnVuY3Rpb24gaW5oZXJpdHMoY3Rvciwgc3VwZXJDdG9yKSB7XG4gICAgY3Rvci5zdXBlcl8gPSBzdXBlckN0b3JcbiAgICBjdG9yLnByb3RvdHlwZSA9IE9iamVjdC5jcmVhdGUoc3VwZXJDdG9yLnByb3RvdHlwZSwge1xuICAgICAgY29uc3RydWN0b3I6IHtcbiAgICAgICAgdmFsdWU6IGN0b3IsXG4gICAgICAgIGVudW1lcmFibGU6IGZhbHNlLFxuICAgICAgICB3cml0YWJsZTogdHJ1ZSxcbiAgICAgICAgY29uZmlndXJhYmxlOiB0cnVlXG4gICAgICB9XG4gICAgfSk7XG4gIH07XG59IGVsc2Uge1xuICAvLyBvbGQgc2Nob29sIHNoaW0gZm9yIG9sZCBicm93c2Vyc1xuICBtb2R1bGUuZXhwb3J0cyA9IGZ1bmN0aW9uIGluaGVyaXRzKGN0b3IsIHN1cGVyQ3Rvcikge1xuICAgIGN0b3Iuc3VwZXJfID0gc3VwZXJDdG9yXG4gICAgdmFyIFRlbXBDdG9yID0gZnVuY3Rpb24gKCkge31cbiAgICBUZW1wQ3Rvci5wcm90b3R5cGUgPSBzdXBlckN0b3IucHJvdG90eXBlXG4gICAgY3Rvci5wcm90b3R5cGUgPSBuZXcgVGVtcEN0b3IoKVxuICAgIGN0b3IucHJvdG90eXBlLmNvbnN0cnVjdG9yID0gY3RvclxuICB9XG59XG4iLCIvLyBzaGltIGZvciB1c2luZyBwcm9jZXNzIGluIGJyb3dzZXJcblxudmFyIHByb2Nlc3MgPSBtb2R1bGUuZXhwb3J0cyA9IHt9O1xuXG5wcm9jZXNzLm5leHRUaWNrID0gKGZ1bmN0aW9uICgpIHtcbiAgICB2YXIgY2FuU2V0SW1tZWRpYXRlID0gdHlwZW9mIHdpbmRvdyAhPT0gJ3VuZGVmaW5lZCdcbiAgICAmJiB3aW5kb3cuc2V0SW1tZWRpYXRlO1xuICAgIHZhciBjYW5Qb3N0ID0gdHlwZW9mIHdpbmRvdyAhPT0gJ3VuZGVmaW5lZCdcbiAgICAmJiB3aW5kb3cucG9zdE1lc3NhZ2UgJiYgd2luZG93LmFkZEV2ZW50TGlzdGVuZXJcbiAgICA7XG5cbiAgICBpZiAoY2FuU2V0SW1tZWRpYXRlKSB7XG4gICAgICAgIHJldHVybiBmdW5jdGlvbiAoZikgeyByZXR1cm4gd2luZG93LnNldEltbWVkaWF0ZShmKSB9O1xuICAgIH1cblxuICAgIGlmIChjYW5Qb3N0KSB7XG4gICAgICAgIHZhciBxdWV1ZSA9IFtdO1xuICAgICAgICB3aW5kb3cuYWRkRXZlbnRMaXN0ZW5lcignbWVzc2FnZScsIGZ1bmN0aW9uIChldikge1xuICAgICAgICAgICAgdmFyIHNvdXJjZSA9IGV2LnNvdXJjZTtcbiAgICAgICAgICAgIGlmICgoc291cmNlID09PSB3aW5kb3cgfHwgc291cmNlID09PSBudWxsKSAmJiBldi5kYXRhID09PSAncHJvY2Vzcy10aWNrJykge1xuICAgICAgICAgICAgICAgIGV2LnN0b3BQcm9wYWdhdGlvbigpO1xuICAgICAgICAgICAgICAgIGlmIChxdWV1ZS5sZW5ndGggPiAwKSB7XG4gICAgICAgICAgICAgICAgICAgIHZhciBmbiA9IHF1ZXVlLnNoaWZ0KCk7XG4gICAgICAgICAgICAgICAgICAgIGZuKCk7XG4gICAgICAgICAgICAgICAgfVxuICAgICAgICAgICAgfVxuICAgICAgICB9LCB0cnVlKTtcblxuICAgICAgICByZXR1cm4gZnVuY3Rpb24gbmV4dFRpY2soZm4pIHtcbiAgICAgICAgICAgIHF1ZXVlLnB1c2goZm4pO1xuICAgICAgICAgICAgd2luZG93LnBvc3RNZXNzYWdlKCdwcm9jZXNzLXRpY2snLCAnKicpO1xuICAgICAgICB9O1xuICAgIH1cblxuICAgIHJldHVybiBmdW5jdGlvbiBuZXh0VGljayhmbikge1xuICAgICAgICBzZXRUaW1lb3V0KGZuLCAwKTtcbiAgICB9O1xufSkoKTtcblxucHJvY2Vzcy50aXRsZSA9ICdicm93c2VyJztcbnByb2Nlc3MuYnJvd3NlciA9IHRydWU7XG5wcm9jZXNzLmVudiA9IHt9O1xucHJvY2Vzcy5hcmd2ID0gW107XG5cbmZ1bmN0aW9uIG5vb3AoKSB7fVxuXG5wcm9jZXNzLm9uID0gbm9vcDtcbnByb2Nlc3MuYWRkTGlzdGVuZXIgPSBub29wO1xucHJvY2Vzcy5vbmNlID0gbm9vcDtcbnByb2Nlc3Mub2ZmID0gbm9vcDtcbnByb2Nlc3MucmVtb3ZlTGlzdGVuZXIgPSBub29wO1xucHJvY2Vzcy5yZW1vdmVBbGxMaXN0ZW5lcnMgPSBub29wO1xucHJvY2Vzcy5lbWl0ID0gbm9vcDtcblxucHJvY2Vzcy5iaW5kaW5nID0gZnVuY3Rpb24gKG5hbWUpIHtcbiAgICB0aHJvdyBuZXcgRXJyb3IoJ3Byb2Nlc3MuYmluZGluZyBpcyBub3Qgc3VwcG9ydGVkJyk7XG59XG5cbi8vIFRPRE8oc2h0eWxtYW4pXG5wcm9jZXNzLmN3ZCA9IGZ1bmN0aW9uICgpIHsgcmV0dXJuICcvJyB9O1xucHJvY2Vzcy5jaGRpciA9IGZ1bmN0aW9uIChkaXIpIHtcbiAgICB0aHJvdyBuZXcgRXJyb3IoJ3Byb2Nlc3MuY2hkaXIgaXMgbm90IHN1cHBvcnRlZCcpO1xufTtcbiIsIm1vZHVsZS5leHBvcnRzID0gZnVuY3Rpb24gaXNCdWZmZXIoYXJnKSB7XG4gIHJldHVybiBhcmcgJiYgdHlwZW9mIGFyZyA9PT0gJ29iamVjdCdcbiAgICAmJiB0eXBlb2YgYXJnLmNvcHkgPT09ICdmdW5jdGlvbidcbiAgICAmJiB0eXBlb2YgYXJnLmZpbGwgPT09ICdmdW5jdGlvbidcbiAgICAmJiB0eXBlb2YgYXJnLnJlYWRVSW50OCA9PT0gJ2Z1bmN0aW9uJztcbn0iLCIoZnVuY3Rpb24gKHByb2Nlc3MsZ2xvYmFsKXtcbi8vIENvcHlyaWdodCBKb3llbnQsIEluYy4gYW5kIG90aGVyIE5vZGUgY29udHJpYnV0b3JzLlxuLy9cbi8vIFBlcm1pc3Npb24gaXMgaGVyZWJ5IGdyYW50ZWQsIGZyZWUgb2YgY2hhcmdlLCB0byBhbnkgcGVyc29uIG9idGFpbmluZyBhXG4vLyBjb3B5IG9mIHRoaXMgc29mdHdhcmUgYW5kIGFzc29jaWF0ZWQgZG9jdW1lbnRhdGlvbiBmaWxlcyAodGhlXG4vLyBcIlNvZnR3YXJlXCIpLCB0byBkZWFsIGluIHRoZSBTb2Z0d2FyZSB3aXRob3V0IHJlc3RyaWN0aW9uLCBpbmNsdWRpbmdcbi8vIHdpdGhvdXQgbGltaXRhdGlvbiB0aGUgcmlnaHRzIHRvIHVzZSwgY29weSwgbW9kaWZ5LCBtZXJnZSwgcHVibGlzaCxcbi8vIGRpc3RyaWJ1dGUsIHN1YmxpY2Vuc2UsIGFuZC9vciBzZWxsIGNvcGllcyBvZiB0aGUgU29mdHdhcmUsIGFuZCB0byBwZXJtaXRcbi8vIHBlcnNvbnMgdG8gd2hvbSB0aGUgU29mdHdhcmUgaXMgZnVybmlzaGVkIHRvIGRvIHNvLCBzdWJqZWN0IHRvIHRoZVxuLy8gZm9sbG93aW5nIGNvbmRpdGlvbnM6XG4vL1xuLy8gVGhlIGFib3ZlIGNvcHlyaWdodCBub3RpY2UgYW5kIHRoaXMgcGVybWlzc2lvbiBub3RpY2Ugc2hhbGwgYmUgaW5jbHVkZWRcbi8vIGluIGFsbCBjb3BpZXMgb3Igc3Vic3RhbnRpYWwgcG9ydGlvbnMgb2YgdGhlIFNvZnR3YXJlLlxuLy9cbi8vIFRIRSBTT0ZUV0FSRSBJUyBQUk9WSURFRCBcIkFTIElTXCIsIFdJVEhPVVQgV0FSUkFOVFkgT0YgQU5ZIEtJTkQsIEVYUFJFU1Ncbi8vIE9SIElNUExJRUQsIElOQ0xVRElORyBCVVQgTk9UIExJTUlURUQgVE8gVEhFIFdBUlJBTlRJRVMgT0Zcbi8vIE1FUkNIQU5UQUJJTElUWSwgRklUTkVTUyBGT1IgQSBQQVJUSUNVTEFSIFBVUlBPU0UgQU5EIE5PTklORlJJTkdFTUVOVC4gSU5cbi8vIE5PIEVWRU5UIFNIQUxMIFRIRSBBVVRIT1JTIE9SIENPUFlSSUdIVCBIT0xERVJTIEJFIExJQUJMRSBGT1IgQU5ZIENMQUlNLFxuLy8gREFNQUdFUyBPUiBPVEhFUiBMSUFCSUxJVFksIFdIRVRIRVIgSU4gQU4gQUNUSU9OIE9GIENPTlRSQUNULCBUT1JUIE9SXG4vLyBPVEhFUldJU0UsIEFSSVNJTkcgRlJPTSwgT1VUIE9GIE9SIElOIENPTk5FQ1RJT04gV0lUSCBUSEUgU09GVFdBUkUgT1IgVEhFXG4vLyBVU0UgT1IgT1RIRVIgREVBTElOR1MgSU4gVEhFIFNPRlRXQVJFLlxuXG52YXIgZm9ybWF0UmVnRXhwID0gLyVbc2RqJV0vZztcbmV4cG9ydHMuZm9ybWF0ID0gZnVuY3Rpb24oZikge1xuICBpZiAoIWlzU3RyaW5nKGYpKSB7XG4gICAgdmFyIG9iamVjdHMgPSBbXTtcbiAgICBmb3IgKHZhciBpID0gMDsgaSA8IGFyZ3VtZW50cy5sZW5ndGg7IGkrKykge1xuICAgICAgb2JqZWN0cy5wdXNoKGluc3BlY3QoYXJndW1lbnRzW2ldKSk7XG4gICAgfVxuICAgIHJldHVybiBvYmplY3RzLmpvaW4oJyAnKTtcbiAgfVxuXG4gIHZhciBpID0gMTtcbiAgdmFyIGFyZ3MgPSBhcmd1bWVudHM7XG4gIHZhciBsZW4gPSBhcmdzLmxlbmd0aDtcbiAgdmFyIHN0ciA9IFN0cmluZyhmKS5yZXBsYWNlKGZvcm1hdFJlZ0V4cCwgZnVuY3Rpb24oeCkge1xuICAgIGlmICh4ID09PSAnJSUnKSByZXR1cm4gJyUnO1xuICAgIGlmIChpID49IGxlbikgcmV0dXJuIHg7XG4gICAgc3dpdGNoICh4KSB7XG4gICAgICBjYXNlICclcyc6IHJldHVybiBTdHJpbmcoYXJnc1tpKytdKTtcbiAgICAgIGNhc2UgJyVkJzogcmV0dXJuIE51bWJlcihhcmdzW2krK10pO1xuICAgICAgY2FzZSAnJWonOlxuICAgICAgICB0cnkge1xuICAgICAgICAgIHJldHVybiBKU09OLnN0cmluZ2lmeShhcmdzW2krK10pO1xuICAgICAgICB9IGNhdGNoIChfKSB7XG4gICAgICAgICAgcmV0dXJuICdbQ2lyY3VsYXJdJztcbiAgICAgICAgfVxuICAgICAgZGVmYXVsdDpcbiAgICAgICAgcmV0dXJuIHg7XG4gICAgfVxuICB9KTtcbiAgZm9yICh2YXIgeCA9IGFyZ3NbaV07IGkgPCBsZW47IHggPSBhcmdzWysraV0pIHtcbiAgICBpZiAoaXNOdWxsKHgpIHx8ICFpc09iamVjdCh4KSkge1xuICAgICAgc3RyICs9ICcgJyArIHg7XG4gICAgfSBlbHNlIHtcbiAgICAgIHN0ciArPSAnICcgKyBpbnNwZWN0KHgpO1xuICAgIH1cbiAgfVxuICByZXR1cm4gc3RyO1xufTtcblxuXG4vLyBNYXJrIHRoYXQgYSBtZXRob2Qgc2hvdWxkIG5vdCBiZSB1c2VkLlxuLy8gUmV0dXJucyBhIG1vZGlmaWVkIGZ1bmN0aW9uIHdoaWNoIHdhcm5zIG9uY2UgYnkgZGVmYXVsdC5cbi8vIElmIC0tbm8tZGVwcmVjYXRpb24gaXMgc2V0LCB0aGVuIGl0IGlzIGEgbm8tb3AuXG5leHBvcnRzLmRlcHJlY2F0ZSA9IGZ1bmN0aW9uKGZuLCBtc2cpIHtcbiAgLy8gQWxsb3cgZm9yIGRlcHJlY2F0aW5nIHRoaW5ncyBpbiB0aGUgcHJvY2VzcyBvZiBzdGFydGluZyB1cC5cbiAgaWYgKGlzVW5kZWZpbmVkKGdsb2JhbC5wcm9jZXNzKSkge1xuICAgIHJldHVybiBmdW5jdGlvbigpIHtcbiAgICAgIHJldHVybiBleHBvcnRzLmRlcHJlY2F0ZShmbiwgbXNnKS5hcHBseSh0aGlzLCBhcmd1bWVudHMpO1xuICAgIH07XG4gIH1cblxuICBpZiAocHJvY2Vzcy5ub0RlcHJlY2F0aW9uID09PSB0cnVlKSB7XG4gICAgcmV0dXJuIGZuO1xuICB9XG5cbiAgdmFyIHdhcm5lZCA9IGZhbHNlO1xuICBmdW5jdGlvbiBkZXByZWNhdGVkKCkge1xuICAgIGlmICghd2FybmVkKSB7XG4gICAgICBpZiAocHJvY2Vzcy50aHJvd0RlcHJlY2F0aW9uKSB7XG4gICAgICAgIHRocm93IG5ldyBFcnJvcihtc2cpO1xuICAgICAgfSBlbHNlIGlmIChwcm9jZXNzLnRyYWNlRGVwcmVjYXRpb24pIHtcbiAgICAgICAgY29uc29sZS50cmFjZShtc2cpO1xuICAgICAgfSBlbHNlIHtcbiAgICAgICAgY29uc29sZS5lcnJvcihtc2cpO1xuICAgICAgfVxuICAgICAgd2FybmVkID0gdHJ1ZTtcbiAgICB9XG4gICAgcmV0dXJuIGZuLmFwcGx5KHRoaXMsIGFyZ3VtZW50cyk7XG4gIH1cblxuICByZXR1cm4gZGVwcmVjYXRlZDtcbn07XG5cblxudmFyIGRlYnVncyA9IHt9O1xudmFyIGRlYnVnRW52aXJvbjtcbmV4cG9ydHMuZGVidWdsb2cgPSBmdW5jdGlvbihzZXQpIHtcbiAgaWYgKGlzVW5kZWZpbmVkKGRlYnVnRW52aXJvbikpXG4gICAgZGVidWdFbnZpcm9uID0gcHJvY2Vzcy5lbnYuTk9ERV9ERUJVRyB8fCAnJztcbiAgc2V0ID0gc2V0LnRvVXBwZXJDYXNlKCk7XG4gIGlmICghZGVidWdzW3NldF0pIHtcbiAgICBpZiAobmV3IFJlZ0V4cCgnXFxcXGInICsgc2V0ICsgJ1xcXFxiJywgJ2knKS50ZXN0KGRlYnVnRW52aXJvbikpIHtcbiAgICAgIHZhciBwaWQgPSBwcm9jZXNzLnBpZDtcbiAgICAgIGRlYnVnc1tzZXRdID0gZnVuY3Rpb24oKSB7XG4gICAgICAgIHZhciBtc2cgPSBleHBvcnRzLmZvcm1hdC5hcHBseShleHBvcnRzLCBhcmd1bWVudHMpO1xuICAgICAgICBjb25zb2xlLmVycm9yKCclcyAlZDogJXMnLCBzZXQsIHBpZCwgbXNnKTtcbiAgICAgIH07XG4gICAgfSBlbHNlIHtcbiAgICAgIGRlYnVnc1tzZXRdID0gZnVuY3Rpb24oKSB7fTtcbiAgICB9XG4gIH1cbiAgcmV0dXJuIGRlYnVnc1tzZXRdO1xufTtcblxuXG4vKipcbiAqIEVjaG9zIHRoZSB2YWx1ZSBvZiBhIHZhbHVlLiBUcnlzIHRvIHByaW50IHRoZSB2YWx1ZSBvdXRcbiAqIGluIHRoZSBiZXN0IHdheSBwb3NzaWJsZSBnaXZlbiB0aGUgZGlmZmVyZW50IHR5cGVzLlxuICpcbiAqIEBwYXJhbSB7T2JqZWN0fSBvYmogVGhlIG9iamVjdCB0byBwcmludCBvdXQuXG4gKiBAcGFyYW0ge09iamVjdH0gb3B0cyBPcHRpb25hbCBvcHRpb25zIG9iamVjdCB0aGF0IGFsdGVycyB0aGUgb3V0cHV0LlxuICovXG4vKiBsZWdhY3k6IG9iaiwgc2hvd0hpZGRlbiwgZGVwdGgsIGNvbG9ycyovXG5mdW5jdGlvbiBpbnNwZWN0KG9iaiwgb3B0cykge1xuICAvLyBkZWZhdWx0IG9wdGlvbnNcbiAgdmFyIGN0eCA9IHtcbiAgICBzZWVuOiBbXSxcbiAgICBzdHlsaXplOiBzdHlsaXplTm9Db2xvclxuICB9O1xuICAvLyBsZWdhY3kuLi5cbiAgaWYgKGFyZ3VtZW50cy5sZW5ndGggPj0gMykgY3R4LmRlcHRoID0gYXJndW1lbnRzWzJdO1xuICBpZiAoYXJndW1lbnRzLmxlbmd0aCA+PSA0KSBjdHguY29sb3JzID0gYXJndW1lbnRzWzNdO1xuICBpZiAoaXNCb29sZWFuKG9wdHMpKSB7XG4gICAgLy8gbGVnYWN5Li4uXG4gICAgY3R4LnNob3dIaWRkZW4gPSBvcHRzO1xuICB9IGVsc2UgaWYgKG9wdHMpIHtcbiAgICAvLyBnb3QgYW4gXCJvcHRpb25zXCIgb2JqZWN0XG4gICAgZXhwb3J0cy5fZXh0ZW5kKGN0eCwgb3B0cyk7XG4gIH1cbiAgLy8gc2V0IGRlZmF1bHQgb3B0aW9uc1xuICBpZiAoaXNVbmRlZmluZWQoY3R4LnNob3dIaWRkZW4pKSBjdHguc2hvd0hpZGRlbiA9IGZhbHNlO1xuICBpZiAoaXNVbmRlZmluZWQoY3R4LmRlcHRoKSkgY3R4LmRlcHRoID0gMjtcbiAgaWYgKGlzVW5kZWZpbmVkKGN0eC5jb2xvcnMpKSBjdHguY29sb3JzID0gZmFsc2U7XG4gIGlmIChpc1VuZGVmaW5lZChjdHguY3VzdG9tSW5zcGVjdCkpIGN0eC5jdXN0b21JbnNwZWN0ID0gdHJ1ZTtcbiAgaWYgKGN0eC5jb2xvcnMpIGN0eC5zdHlsaXplID0gc3R5bGl6ZVdpdGhDb2xvcjtcbiAgcmV0dXJuIGZvcm1hdFZhbHVlKGN0eCwgb2JqLCBjdHguZGVwdGgpO1xufVxuZXhwb3J0cy5pbnNwZWN0ID0gaW5zcGVjdDtcblxuXG4vLyBodHRwOi8vZW4ud2lraXBlZGlhLm9yZy93aWtpL0FOU0lfZXNjYXBlX2NvZGUjZ3JhcGhpY3Ncbmluc3BlY3QuY29sb3JzID0ge1xuICAnYm9sZCcgOiBbMSwgMjJdLFxuICAnaXRhbGljJyA6IFszLCAyM10sXG4gICd1bmRlcmxpbmUnIDogWzQsIDI0XSxcbiAgJ2ludmVyc2UnIDogWzcsIDI3XSxcbiAgJ3doaXRlJyA6IFszNywgMzldLFxuICAnZ3JleScgOiBbOTAsIDM5XSxcbiAgJ2JsYWNrJyA6IFszMCwgMzldLFxuICAnYmx1ZScgOiBbMzQsIDM5XSxcbiAgJ2N5YW4nIDogWzM2LCAzOV0sXG4gICdncmVlbicgOiBbMzIsIDM5XSxcbiAgJ21hZ2VudGEnIDogWzM1LCAzOV0sXG4gICdyZWQnIDogWzMxLCAzOV0sXG4gICd5ZWxsb3cnIDogWzMzLCAzOV1cbn07XG5cbi8vIERvbid0IHVzZSAnYmx1ZScgbm90IHZpc2libGUgb24gY21kLmV4ZVxuaW5zcGVjdC5zdHlsZXMgPSB7XG4gICdzcGVjaWFsJzogJ2N5YW4nLFxuICAnbnVtYmVyJzogJ3llbGxvdycsXG4gICdib29sZWFuJzogJ3llbGxvdycsXG4gICd1bmRlZmluZWQnOiAnZ3JleScsXG4gICdudWxsJzogJ2JvbGQnLFxuICAnc3RyaW5nJzogJ2dyZWVuJyxcbiAgJ2RhdGUnOiAnbWFnZW50YScsXG4gIC8vIFwibmFtZVwiOiBpbnRlbnRpb25hbGx5IG5vdCBzdHlsaW5nXG4gICdyZWdleHAnOiAncmVkJ1xufTtcblxuXG5mdW5jdGlvbiBzdHlsaXplV2l0aENvbG9yKHN0ciwgc3R5bGVUeXBlKSB7XG4gIHZhciBzdHlsZSA9IGluc3BlY3Quc3R5bGVzW3N0eWxlVHlwZV07XG5cbiAgaWYgKHN0eWxlKSB7XG4gICAgcmV0dXJuICdcXHUwMDFiWycgKyBpbnNwZWN0LmNvbG9yc1tzdHlsZV1bMF0gKyAnbScgKyBzdHIgK1xuICAgICAgICAgICAnXFx1MDAxYlsnICsgaW5zcGVjdC5jb2xvcnNbc3R5bGVdWzFdICsgJ20nO1xuICB9IGVsc2Uge1xuICAgIHJldHVybiBzdHI7XG4gIH1cbn1cblxuXG5mdW5jdGlvbiBzdHlsaXplTm9Db2xvcihzdHIsIHN0eWxlVHlwZSkge1xuICByZXR1cm4gc3RyO1xufVxuXG5cbmZ1bmN0aW9uIGFycmF5VG9IYXNoKGFycmF5KSB7XG4gIHZhciBoYXNoID0ge307XG5cbiAgYXJyYXkuZm9yRWFjaChmdW5jdGlvbih2YWwsIGlkeCkge1xuICAgIGhhc2hbdmFsXSA9IHRydWU7XG4gIH0pO1xuXG4gIHJldHVybiBoYXNoO1xufVxuXG5cbmZ1bmN0aW9uIGZvcm1hdFZhbHVlKGN0eCwgdmFsdWUsIHJlY3Vyc2VUaW1lcykge1xuICAvLyBQcm92aWRlIGEgaG9vayBmb3IgdXNlci1zcGVjaWZpZWQgaW5zcGVjdCBmdW5jdGlvbnMuXG4gIC8vIENoZWNrIHRoYXQgdmFsdWUgaXMgYW4gb2JqZWN0IHdpdGggYW4gaW5zcGVjdCBmdW5jdGlvbiBvbiBpdFxuICBpZiAoY3R4LmN1c3RvbUluc3BlY3QgJiZcbiAgICAgIHZhbHVlICYmXG4gICAgICBpc0Z1bmN0aW9uKHZhbHVlLmluc3BlY3QpICYmXG4gICAgICAvLyBGaWx0ZXIgb3V0IHRoZSB1dGlsIG1vZHVsZSwgaXQncyBpbnNwZWN0IGZ1bmN0aW9uIGlzIHNwZWNpYWxcbiAgICAgIHZhbHVlLmluc3BlY3QgIT09IGV4cG9ydHMuaW5zcGVjdCAmJlxuICAgICAgLy8gQWxzbyBmaWx0ZXIgb3V0IGFueSBwcm90b3R5cGUgb2JqZWN0cyB1c2luZyB0aGUgY2lyY3VsYXIgY2hlY2suXG4gICAgICAhKHZhbHVlLmNvbnN0cnVjdG9yICYmIHZhbHVlLmNvbnN0cnVjdG9yLnByb3RvdHlwZSA9PT0gdmFsdWUpKSB7XG4gICAgdmFyIHJldCA9IHZhbHVlLmluc3BlY3QocmVjdXJzZVRpbWVzLCBjdHgpO1xuICAgIGlmICghaXNTdHJpbmcocmV0KSkge1xuICAgICAgcmV0ID0gZm9ybWF0VmFsdWUoY3R4LCByZXQsIHJlY3Vyc2VUaW1lcyk7XG4gICAgfVxuICAgIHJldHVybiByZXQ7XG4gIH1cblxuICAvLyBQcmltaXRpdmUgdHlwZXMgY2Fubm90IGhhdmUgcHJvcGVydGllc1xuICB2YXIgcHJpbWl0aXZlID0gZm9ybWF0UHJpbWl0aXZlKGN0eCwgdmFsdWUpO1xuICBpZiAocHJpbWl0aXZlKSB7XG4gICAgcmV0dXJuIHByaW1pdGl2ZTtcbiAgfVxuXG4gIC8vIExvb2sgdXAgdGhlIGtleXMgb2YgdGhlIG9iamVjdC5cbiAgdmFyIGtleXMgPSBPYmplY3Qua2V5cyh2YWx1ZSk7XG4gIHZhciB2aXNpYmxlS2V5cyA9IGFycmF5VG9IYXNoKGtleXMpO1xuXG4gIGlmIChjdHguc2hvd0hpZGRlbikge1xuICAgIGtleXMgPSBPYmplY3QuZ2V0T3duUHJvcGVydHlOYW1lcyh2YWx1ZSk7XG4gIH1cblxuICAvLyBJRSBkb2Vzbid0IG1ha2UgZXJyb3IgZmllbGRzIG5vbi1lbnVtZXJhYmxlXG4gIC8vIGh0dHA6Ly9tc2RuLm1pY3Jvc29mdC5jb20vZW4tdXMvbGlicmFyeS9pZS9kd3c1MnNidCh2PXZzLjk0KS5hc3B4XG4gIGlmIChpc0Vycm9yKHZhbHVlKVxuICAgICAgJiYgKGtleXMuaW5kZXhPZignbWVzc2FnZScpID49IDAgfHwga2V5cy5pbmRleE9mKCdkZXNjcmlwdGlvbicpID49IDApKSB7XG4gICAgcmV0dXJuIGZvcm1hdEVycm9yKHZhbHVlKTtcbiAgfVxuXG4gIC8vIFNvbWUgdHlwZSBvZiBvYmplY3Qgd2l0aG91dCBwcm9wZXJ0aWVzIGNhbiBiZSBzaG9ydGN1dHRlZC5cbiAgaWYgKGtleXMubGVuZ3RoID09PSAwKSB7XG4gICAgaWYgKGlzRnVuY3Rpb24odmFsdWUpKSB7XG4gICAgICB2YXIgbmFtZSA9IHZhbHVlLm5hbWUgPyAnOiAnICsgdmFsdWUubmFtZSA6ICcnO1xuICAgICAgcmV0dXJuIGN0eC5zdHlsaXplKCdbRnVuY3Rpb24nICsgbmFtZSArICddJywgJ3NwZWNpYWwnKTtcbiAgICB9XG4gICAgaWYgKGlzUmVnRXhwKHZhbHVlKSkge1xuICAgICAgcmV0dXJuIGN0eC5zdHlsaXplKFJlZ0V4cC5wcm90b3R5cGUudG9TdHJpbmcuY2FsbCh2YWx1ZSksICdyZWdleHAnKTtcbiAgICB9XG4gICAgaWYgKGlzRGF0ZSh2YWx1ZSkpIHtcbiAgICAgIHJldHVybiBjdHguc3R5bGl6ZShEYXRlLnByb3RvdHlwZS50b1N0cmluZy5jYWxsKHZhbHVlKSwgJ2RhdGUnKTtcbiAgICB9XG4gICAgaWYgKGlzRXJyb3IodmFsdWUpKSB7XG4gICAgICByZXR1cm4gZm9ybWF0RXJyb3IodmFsdWUpO1xuICAgIH1cbiAgfVxuXG4gIHZhciBiYXNlID0gJycsIGFycmF5ID0gZmFsc2UsIGJyYWNlcyA9IFsneycsICd9J107XG5cbiAgLy8gTWFrZSBBcnJheSBzYXkgdGhhdCB0aGV5IGFyZSBBcnJheVxuICBpZiAoaXNBcnJheSh2YWx1ZSkpIHtcbiAgICBhcnJheSA9IHRydWU7XG4gICAgYnJhY2VzID0gWydbJywgJ10nXTtcbiAgfVxuXG4gIC8vIE1ha2UgZnVuY3Rpb25zIHNheSB0aGF0IHRoZXkgYXJlIGZ1bmN0aW9uc1xuICBpZiAoaXNGdW5jdGlvbih2YWx1ZSkpIHtcbiAgICB2YXIgbiA9IHZhbHVlLm5hbWUgPyAnOiAnICsgdmFsdWUubmFtZSA6ICcnO1xuICAgIGJhc2UgPSAnIFtGdW5jdGlvbicgKyBuICsgJ10nO1xuICB9XG5cbiAgLy8gTWFrZSBSZWdFeHBzIHNheSB0aGF0IHRoZXkgYXJlIFJlZ0V4cHNcbiAgaWYgKGlzUmVnRXhwKHZhbHVlKSkge1xuICAgIGJhc2UgPSAnICcgKyBSZWdFeHAucHJvdG90eXBlLnRvU3RyaW5nLmNhbGwodmFsdWUpO1xuICB9XG5cbiAgLy8gTWFrZSBkYXRlcyB3aXRoIHByb3BlcnRpZXMgZmlyc3Qgc2F5IHRoZSBkYXRlXG4gIGlmIChpc0RhdGUodmFsdWUpKSB7XG4gICAgYmFzZSA9ICcgJyArIERhdGUucHJvdG90eXBlLnRvVVRDU3RyaW5nLmNhbGwodmFsdWUpO1xuICB9XG5cbiAgLy8gTWFrZSBlcnJvciB3aXRoIG1lc3NhZ2UgZmlyc3Qgc2F5IHRoZSBlcnJvclxuICBpZiAoaXNFcnJvcih2YWx1ZSkpIHtcbiAgICBiYXNlID0gJyAnICsgZm9ybWF0RXJyb3IodmFsdWUpO1xuICB9XG5cbiAgaWYgKGtleXMubGVuZ3RoID09PSAwICYmICghYXJyYXkgfHwgdmFsdWUubGVuZ3RoID09IDApKSB7XG4gICAgcmV0dXJuIGJyYWNlc1swXSArIGJhc2UgKyBicmFjZXNbMV07XG4gIH1cblxuICBpZiAocmVjdXJzZVRpbWVzIDwgMCkge1xuICAgIGlmIChpc1JlZ0V4cCh2YWx1ZSkpIHtcbiAgICAgIHJldHVybiBjdHguc3R5bGl6ZShSZWdFeHAucHJvdG90eXBlLnRvU3RyaW5nLmNhbGwodmFsdWUpLCAncmVnZXhwJyk7XG4gICAgfSBlbHNlIHtcbiAgICAgIHJldHVybiBjdHguc3R5bGl6ZSgnW09iamVjdF0nLCAnc3BlY2lhbCcpO1xuICAgIH1cbiAgfVxuXG4gIGN0eC5zZWVuLnB1c2godmFsdWUpO1xuXG4gIHZhciBvdXRwdXQ7XG4gIGlmIChhcnJheSkge1xuICAgIG91dHB1dCA9IGZvcm1hdEFycmF5KGN0eCwgdmFsdWUsIHJlY3Vyc2VUaW1lcywgdmlzaWJsZUtleXMsIGtleXMpO1xuICB9IGVsc2Uge1xuICAgIG91dHB1dCA9IGtleXMubWFwKGZ1bmN0aW9uKGtleSkge1xuICAgICAgcmV0dXJuIGZvcm1hdFByb3BlcnR5KGN0eCwgdmFsdWUsIHJlY3Vyc2VUaW1lcywgdmlzaWJsZUtleXMsIGtleSwgYXJyYXkpO1xuICAgIH0pO1xuICB9XG5cbiAgY3R4LnNlZW4ucG9wKCk7XG5cbiAgcmV0dXJuIHJlZHVjZVRvU2luZ2xlU3RyaW5nKG91dHB1dCwgYmFzZSwgYnJhY2VzKTtcbn1cblxuXG5mdW5jdGlvbiBmb3JtYXRQcmltaXRpdmUoY3R4LCB2YWx1ZSkge1xuICBpZiAoaXNVbmRlZmluZWQodmFsdWUpKVxuICAgIHJldHVybiBjdHguc3R5bGl6ZSgndW5kZWZpbmVkJywgJ3VuZGVmaW5lZCcpO1xuICBpZiAoaXNTdHJpbmcodmFsdWUpKSB7XG4gICAgdmFyIHNpbXBsZSA9ICdcXCcnICsgSlNPTi5zdHJpbmdpZnkodmFsdWUpLnJlcGxhY2UoL15cInxcIiQvZywgJycpXG4gICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAucmVwbGFjZSgvJy9nLCBcIlxcXFwnXCIpXG4gICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAucmVwbGFjZSgvXFxcXFwiL2csICdcIicpICsgJ1xcJyc7XG4gICAgcmV0dXJuIGN0eC5zdHlsaXplKHNpbXBsZSwgJ3N0cmluZycpO1xuICB9XG4gIGlmIChpc051bWJlcih2YWx1ZSkpXG4gICAgcmV0dXJuIGN0eC5zdHlsaXplKCcnICsgdmFsdWUsICdudW1iZXInKTtcbiAgaWYgKGlzQm9vbGVhbih2YWx1ZSkpXG4gICAgcmV0dXJuIGN0eC5zdHlsaXplKCcnICsgdmFsdWUsICdib29sZWFuJyk7XG4gIC8vIEZvciBzb21lIHJlYXNvbiB0eXBlb2YgbnVsbCBpcyBcIm9iamVjdFwiLCBzbyBzcGVjaWFsIGNhc2UgaGVyZS5cbiAgaWYgKGlzTnVsbCh2YWx1ZSkpXG4gICAgcmV0dXJuIGN0eC5zdHlsaXplKCdudWxsJywgJ251bGwnKTtcbn1cblxuXG5mdW5jdGlvbiBmb3JtYXRFcnJvcih2YWx1ZSkge1xuICByZXR1cm4gJ1snICsgRXJyb3IucHJvdG90eXBlLnRvU3RyaW5nLmNhbGwodmFsdWUpICsgJ10nO1xufVxuXG5cbmZ1bmN0aW9uIGZvcm1hdEFycmF5KGN0eCwgdmFsdWUsIHJlY3Vyc2VUaW1lcywgdmlzaWJsZUtleXMsIGtleXMpIHtcbiAgdmFyIG91dHB1dCA9IFtdO1xuICBmb3IgKHZhciBpID0gMCwgbCA9IHZhbHVlLmxlbmd0aDsgaSA8IGw7ICsraSkge1xuICAgIGlmIChoYXNPd25Qcm9wZXJ0eSh2YWx1ZSwgU3RyaW5nKGkpKSkge1xuICAgICAgb3V0cHV0LnB1c2goZm9ybWF0UHJvcGVydHkoY3R4LCB2YWx1ZSwgcmVjdXJzZVRpbWVzLCB2aXNpYmxlS2V5cyxcbiAgICAgICAgICBTdHJpbmcoaSksIHRydWUpKTtcbiAgICB9IGVsc2Uge1xuICAgICAgb3V0cHV0LnB1c2goJycpO1xuICAgIH1cbiAgfVxuICBrZXlzLmZvckVhY2goZnVuY3Rpb24oa2V5KSB7XG4gICAgaWYgKCFrZXkubWF0Y2goL15cXGQrJC8pKSB7XG4gICAgICBvdXRwdXQucHVzaChmb3JtYXRQcm9wZXJ0eShjdHgsIHZhbHVlLCByZWN1cnNlVGltZXMsIHZpc2libGVLZXlzLFxuICAgICAgICAgIGtleSwgdHJ1ZSkpO1xuICAgIH1cbiAgfSk7XG4gIHJldHVybiBvdXRwdXQ7XG59XG5cblxuZnVuY3Rpb24gZm9ybWF0UHJvcGVydHkoY3R4LCB2YWx1ZSwgcmVjdXJzZVRpbWVzLCB2aXNpYmxlS2V5cywga2V5LCBhcnJheSkge1xuICB2YXIgbmFtZSwgc3RyLCBkZXNjO1xuICBkZXNjID0gT2JqZWN0LmdldE93blByb3BlcnR5RGVzY3JpcHRvcih2YWx1ZSwga2V5KSB8fCB7IHZhbHVlOiB2YWx1ZVtrZXldIH07XG4gIGlmIChkZXNjLmdldCkge1xuICAgIGlmIChkZXNjLnNldCkge1xuICAgICAgc3RyID0gY3R4LnN0eWxpemUoJ1tHZXR0ZXIvU2V0dGVyXScsICdzcGVjaWFsJyk7XG4gICAgfSBlbHNlIHtcbiAgICAgIHN0ciA9IGN0eC5zdHlsaXplKCdbR2V0dGVyXScsICdzcGVjaWFsJyk7XG4gICAgfVxuICB9IGVsc2Uge1xuICAgIGlmIChkZXNjLnNldCkge1xuICAgICAgc3RyID0gY3R4LnN0eWxpemUoJ1tTZXR0ZXJdJywgJ3NwZWNpYWwnKTtcbiAgICB9XG4gIH1cbiAgaWYgKCFoYXNPd25Qcm9wZXJ0eSh2aXNpYmxlS2V5cywga2V5KSkge1xuICAgIG5hbWUgPSAnWycgKyBrZXkgKyAnXSc7XG4gIH1cbiAgaWYgKCFzdHIpIHtcbiAgICBpZiAoY3R4LnNlZW4uaW5kZXhPZihkZXNjLnZhbHVlKSA8IDApIHtcbiAgICAgIGlmIChpc051bGwocmVjdXJzZVRpbWVzKSkge1xuICAgICAgICBzdHIgPSBmb3JtYXRWYWx1ZShjdHgsIGRlc2MudmFsdWUsIG51bGwpO1xuICAgICAgfSBlbHNlIHtcbiAgICAgICAgc3RyID0gZm9ybWF0VmFsdWUoY3R4LCBkZXNjLnZhbHVlLCByZWN1cnNlVGltZXMgLSAxKTtcbiAgICAgIH1cbiAgICAgIGlmIChzdHIuaW5kZXhPZignXFxuJykgPiAtMSkge1xuICAgICAgICBpZiAoYXJyYXkpIHtcbiAgICAgICAgICBzdHIgPSBzdHIuc3BsaXQoJ1xcbicpLm1hcChmdW5jdGlvbihsaW5lKSB7XG4gICAgICAgICAgICByZXR1cm4gJyAgJyArIGxpbmU7XG4gICAgICAgICAgfSkuam9pbignXFxuJykuc3Vic3RyKDIpO1xuICAgICAgICB9IGVsc2Uge1xuICAgICAgICAgIHN0ciA9ICdcXG4nICsgc3RyLnNwbGl0KCdcXG4nKS5tYXAoZnVuY3Rpb24obGluZSkge1xuICAgICAgICAgICAgcmV0dXJuICcgICAnICsgbGluZTtcbiAgICAgICAgICB9KS5qb2luKCdcXG4nKTtcbiAgICAgICAgfVxuICAgICAgfVxuICAgIH0gZWxzZSB7XG4gICAgICBzdHIgPSBjdHguc3R5bGl6ZSgnW0NpcmN1bGFyXScsICdzcGVjaWFsJyk7XG4gICAgfVxuICB9XG4gIGlmIChpc1VuZGVmaW5lZChuYW1lKSkge1xuICAgIGlmIChhcnJheSAmJiBrZXkubWF0Y2goL15cXGQrJC8pKSB7XG4gICAgICByZXR1cm4gc3RyO1xuICAgIH1cbiAgICBuYW1lID0gSlNPTi5zdHJpbmdpZnkoJycgKyBrZXkpO1xuICAgIGlmIChuYW1lLm1hdGNoKC9eXCIoW2EtekEtWl9dW2EtekEtWl8wLTldKilcIiQvKSkge1xuICAgICAgbmFtZSA9IG5hbWUuc3Vic3RyKDEsIG5hbWUubGVuZ3RoIC0gMik7XG4gICAgICBuYW1lID0gY3R4LnN0eWxpemUobmFtZSwgJ25hbWUnKTtcbiAgICB9IGVsc2Uge1xuICAgICAgbmFtZSA9IG5hbWUucmVwbGFjZSgvJy9nLCBcIlxcXFwnXCIpXG4gICAgICAgICAgICAgICAgIC5yZXBsYWNlKC9cXFxcXCIvZywgJ1wiJylcbiAgICAgICAgICAgICAgICAgLnJlcGxhY2UoLyheXCJ8XCIkKS9nLCBcIidcIik7XG4gICAgICBuYW1lID0gY3R4LnN0eWxpemUobmFtZSwgJ3N0cmluZycpO1xuICAgIH1cbiAgfVxuXG4gIHJldHVybiBuYW1lICsgJzogJyArIHN0cjtcbn1cblxuXG5mdW5jdGlvbiByZWR1Y2VUb1NpbmdsZVN0cmluZyhvdXRwdXQsIGJhc2UsIGJyYWNlcykge1xuICB2YXIgbnVtTGluZXNFc3QgPSAwO1xuICB2YXIgbGVuZ3RoID0gb3V0cHV0LnJlZHVjZShmdW5jdGlvbihwcmV2LCBjdXIpIHtcbiAgICBudW1MaW5lc0VzdCsrO1xuICAgIGlmIChjdXIuaW5kZXhPZignXFxuJykgPj0gMCkgbnVtTGluZXNFc3QrKztcbiAgICByZXR1cm4gcHJldiArIGN1ci5yZXBsYWNlKC9cXHUwMDFiXFxbXFxkXFxkP20vZywgJycpLmxlbmd0aCArIDE7XG4gIH0sIDApO1xuXG4gIGlmIChsZW5ndGggPiA2MCkge1xuICAgIHJldHVybiBicmFjZXNbMF0gK1xuICAgICAgICAgICAoYmFzZSA9PT0gJycgPyAnJyA6IGJhc2UgKyAnXFxuICcpICtcbiAgICAgICAgICAgJyAnICtcbiAgICAgICAgICAgb3V0cHV0LmpvaW4oJyxcXG4gICcpICtcbiAgICAgICAgICAgJyAnICtcbiAgICAgICAgICAgYnJhY2VzWzFdO1xuICB9XG5cbiAgcmV0dXJuIGJyYWNlc1swXSArIGJhc2UgKyAnICcgKyBvdXRwdXQuam9pbignLCAnKSArICcgJyArIGJyYWNlc1sxXTtcbn1cblxuXG4vLyBOT1RFOiBUaGVzZSB0eXBlIGNoZWNraW5nIGZ1bmN0aW9ucyBpbnRlbnRpb25hbGx5IGRvbid0IHVzZSBgaW5zdGFuY2VvZmBcbi8vIGJlY2F1c2UgaXQgaXMgZnJhZ2lsZSBhbmQgY2FuIGJlIGVhc2lseSBmYWtlZCB3aXRoIGBPYmplY3QuY3JlYXRlKClgLlxuZnVuY3Rpb24gaXNBcnJheShhcikge1xuICByZXR1cm4gQXJyYXkuaXNBcnJheShhcik7XG59XG5leHBvcnRzLmlzQXJyYXkgPSBpc0FycmF5O1xuXG5mdW5jdGlvbiBpc0Jvb2xlYW4oYXJnKSB7XG4gIHJldHVybiB0eXBlb2YgYXJnID09PSAnYm9vbGVhbic7XG59XG5leHBvcnRzLmlzQm9vbGVhbiA9IGlzQm9vbGVhbjtcblxuZnVuY3Rpb24gaXNOdWxsKGFyZykge1xuICByZXR1cm4gYXJnID09PSBudWxsO1xufVxuZXhwb3J0cy5pc051bGwgPSBpc051bGw7XG5cbmZ1bmN0aW9uIGlzTnVsbE9yVW5kZWZpbmVkKGFyZykge1xuICByZXR1cm4gYXJnID09IG51bGw7XG59XG5leHBvcnRzLmlzTnVsbE9yVW5kZWZpbmVkID0gaXNOdWxsT3JVbmRlZmluZWQ7XG5cbmZ1bmN0aW9uIGlzTnVtYmVyKGFyZykge1xuICByZXR1cm4gdHlwZW9mIGFyZyA9PT0gJ251bWJlcic7XG59XG5leHBvcnRzLmlzTnVtYmVyID0gaXNOdW1iZXI7XG5cbmZ1bmN0aW9uIGlzU3RyaW5nKGFyZykge1xuICByZXR1cm4gdHlwZW9mIGFyZyA9PT0gJ3N0cmluZyc7XG59XG5leHBvcnRzLmlzU3RyaW5nID0gaXNTdHJpbmc7XG5cbmZ1bmN0aW9uIGlzU3ltYm9sKGFyZykge1xuICByZXR1cm4gdHlwZW9mIGFyZyA9PT0gJ3N5bWJvbCc7XG59XG5leHBvcnRzLmlzU3ltYm9sID0gaXNTeW1ib2w7XG5cbmZ1bmN0aW9uIGlzVW5kZWZpbmVkKGFyZykge1xuICByZXR1cm4gYXJnID09PSB2b2lkIDA7XG59XG5leHBvcnRzLmlzVW5kZWZpbmVkID0gaXNVbmRlZmluZWQ7XG5cbmZ1bmN0aW9uIGlzUmVnRXhwKHJlKSB7XG4gIHJldHVybiBpc09iamVjdChyZSkgJiYgb2JqZWN0VG9TdHJpbmcocmUpID09PSAnW29iamVjdCBSZWdFeHBdJztcbn1cbmV4cG9ydHMuaXNSZWdFeHAgPSBpc1JlZ0V4cDtcblxuZnVuY3Rpb24gaXNPYmplY3QoYXJnKSB7XG4gIHJldHVybiB0eXBlb2YgYXJnID09PSAnb2JqZWN0JyAmJiBhcmcgIT09IG51bGw7XG59XG5leHBvcnRzLmlzT2JqZWN0ID0gaXNPYmplY3Q7XG5cbmZ1bmN0aW9uIGlzRGF0ZShkKSB7XG4gIHJldHVybiBpc09iamVjdChkKSAmJiBvYmplY3RUb1N0cmluZyhkKSA9PT0gJ1tvYmplY3QgRGF0ZV0nO1xufVxuZXhwb3J0cy5pc0RhdGUgPSBpc0RhdGU7XG5cbmZ1bmN0aW9uIGlzRXJyb3IoZSkge1xuICByZXR1cm4gaXNPYmplY3QoZSkgJiZcbiAgICAgIChvYmplY3RUb1N0cmluZyhlKSA9PT0gJ1tvYmplY3QgRXJyb3JdJyB8fCBlIGluc3RhbmNlb2YgRXJyb3IpO1xufVxuZXhwb3J0cy5pc0Vycm9yID0gaXNFcnJvcjtcblxuZnVuY3Rpb24gaXNGdW5jdGlvbihhcmcpIHtcbiAgcmV0dXJuIHR5cGVvZiBhcmcgPT09ICdmdW5jdGlvbic7XG59XG5leHBvcnRzLmlzRnVuY3Rpb24gPSBpc0Z1bmN0aW9uO1xuXG5mdW5jdGlvbiBpc1ByaW1pdGl2ZShhcmcpIHtcbiAgcmV0dXJuIGFyZyA9PT0gbnVsbCB8fFxuICAgICAgICAgdHlwZW9mIGFyZyA9PT0gJ2Jvb2xlYW4nIHx8XG4gICAgICAgICB0eXBlb2YgYXJnID09PSAnbnVtYmVyJyB8fFxuICAgICAgICAgdHlwZW9mIGFyZyA9PT0gJ3N0cmluZycgfHxcbiAgICAgICAgIHR5cGVvZiBhcmcgPT09ICdzeW1ib2wnIHx8ICAvLyBFUzYgc3ltYm9sXG4gICAgICAgICB0eXBlb2YgYXJnID09PSAndW5kZWZpbmVkJztcbn1cbmV4cG9ydHMuaXNQcmltaXRpdmUgPSBpc1ByaW1pdGl2ZTtcblxuZXhwb3J0cy5pc0J1ZmZlciA9IHJlcXVpcmUoJy4vc3VwcG9ydC9pc0J1ZmZlcicpO1xuXG5mdW5jdGlvbiBvYmplY3RUb1N0cmluZyhvKSB7XG4gIHJldHVybiBPYmplY3QucHJvdG90eXBlLnRvU3RyaW5nLmNhbGwobyk7XG59XG5cblxuZnVuY3Rpb24gcGFkKG4pIHtcbiAgcmV0dXJuIG4gPCAxMCA/ICcwJyArIG4udG9TdHJpbmcoMTApIDogbi50b1N0cmluZygxMCk7XG59XG5cblxudmFyIG1vbnRocyA9IFsnSmFuJywgJ0ZlYicsICdNYXInLCAnQXByJywgJ01heScsICdKdW4nLCAnSnVsJywgJ0F1ZycsICdTZXAnLFxuICAgICAgICAgICAgICAnT2N0JywgJ05vdicsICdEZWMnXTtcblxuLy8gMjYgRmViIDE2OjE5OjM0XG5mdW5jdGlvbiB0aW1lc3RhbXAoKSB7XG4gIHZhciBkID0gbmV3IERhdGUoKTtcbiAgdmFyIHRpbWUgPSBbcGFkKGQuZ2V0SG91cnMoKSksXG4gICAgICAgICAgICAgIHBhZChkLmdldE1pbnV0ZXMoKSksXG4gICAgICAgICAgICAgIHBhZChkLmdldFNlY29uZHMoKSldLmpvaW4oJzonKTtcbiAgcmV0dXJuIFtkLmdldERhdGUoKSwgbW9udGhzW2QuZ2V0TW9udGgoKV0sIHRpbWVdLmpvaW4oJyAnKTtcbn1cblxuXG4vLyBsb2cgaXMganVzdCBhIHRoaW4gd3JhcHBlciB0byBjb25zb2xlLmxvZyB0aGF0IHByZXBlbmRzIGEgdGltZXN0YW1wXG5leHBvcnRzLmxvZyA9IGZ1bmN0aW9uKCkge1xuICBjb25zb2xlLmxvZygnJXMgLSAlcycsIHRpbWVzdGFtcCgpLCBleHBvcnRzLmZvcm1hdC5hcHBseShleHBvcnRzLCBhcmd1bWVudHMpKTtcbn07XG5cblxuLyoqXG4gKiBJbmhlcml0IHRoZSBwcm90b3R5cGUgbWV0aG9kcyBmcm9tIG9uZSBjb25zdHJ1Y3RvciBpbnRvIGFub3RoZXIuXG4gKlxuICogVGhlIEZ1bmN0aW9uLnByb3RvdHlwZS5pbmhlcml0cyBmcm9tIGxhbmcuanMgcmV3cml0dGVuIGFzIGEgc3RhbmRhbG9uZVxuICogZnVuY3Rpb24gKG5vdCBvbiBGdW5jdGlvbi5wcm90b3R5cGUpLiBOT1RFOiBJZiB0aGlzIGZpbGUgaXMgdG8gYmUgbG9hZGVkXG4gKiBkdXJpbmcgYm9vdHN0cmFwcGluZyB0aGlzIGZ1bmN0aW9uIG5lZWRzIHRvIGJlIHJld3JpdHRlbiB1c2luZyBzb21lIG5hdGl2ZVxuICogZnVuY3Rpb25zIGFzIHByb3RvdHlwZSBzZXR1cCB1c2luZyBub3JtYWwgSmF2YVNjcmlwdCBkb2VzIG5vdCB3b3JrIGFzXG4gKiBleHBlY3RlZCBkdXJpbmcgYm9vdHN0cmFwcGluZyAoc2VlIG1pcnJvci5qcyBpbiByMTE0OTAzKS5cbiAqXG4gKiBAcGFyYW0ge2Z1bmN0aW9ufSBjdG9yIENvbnN0cnVjdG9yIGZ1bmN0aW9uIHdoaWNoIG5lZWRzIHRvIGluaGVyaXQgdGhlXG4gKiAgICAgcHJvdG90eXBlLlxuICogQHBhcmFtIHtmdW5jdGlvbn0gc3VwZXJDdG9yIENvbnN0cnVjdG9yIGZ1bmN0aW9uIHRvIGluaGVyaXQgcHJvdG90eXBlIGZyb20uXG4gKi9cbmV4cG9ydHMuaW5oZXJpdHMgPSByZXF1aXJlKCdpbmhlcml0cycpO1xuXG5leHBvcnRzLl9leHRlbmQgPSBmdW5jdGlvbihvcmlnaW4sIGFkZCkge1xuICAvLyBEb24ndCBkbyBhbnl0aGluZyBpZiBhZGQgaXNuJ3QgYW4gb2JqZWN0XG4gIGlmICghYWRkIHx8ICFpc09iamVjdChhZGQpKSByZXR1cm4gb3JpZ2luO1xuXG4gIHZhciBrZXlzID0gT2JqZWN0LmtleXMoYWRkKTtcbiAgdmFyIGkgPSBrZXlzLmxlbmd0aDtcbiAgd2hpbGUgKGktLSkge1xuICAgIG9yaWdpbltrZXlzW2ldXSA9IGFkZFtrZXlzW2ldXTtcbiAgfVxuICByZXR1cm4gb3JpZ2luO1xufTtcblxuZnVuY3Rpb24gaGFzT3duUHJvcGVydHkob2JqLCBwcm9wKSB7XG4gIHJldHVybiBPYmplY3QucHJvdG90eXBlLmhhc093blByb3BlcnR5LmNhbGwob2JqLCBwcm9wKTtcbn1cblxufSkuY2FsbCh0aGlzLHJlcXVpcmUoJ19wcm9jZXNzJyksdHlwZW9mIGdsb2JhbCAhPT0gXCJ1bmRlZmluZWRcIiA/IGdsb2JhbCA6IHR5cGVvZiBzZWxmICE9PSBcInVuZGVmaW5lZFwiID8gc2VsZiA6IHR5cGVvZiB3aW5kb3cgIT09IFwidW5kZWZpbmVkXCIgPyB3aW5kb3cgOiB7fSkiXX0=
