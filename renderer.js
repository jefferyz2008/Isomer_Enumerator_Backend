/**
 * This file is loaded via the <script> tag in the index.html file and will
 * be executed in the renderer process for that window. No Node.js APIs are
 * available in this process because `nodeIntegration` is turned off and
 * `contextIsolation` is turned on. Use the contextBridge API in `preload.js`
 * to expose Node.js functionality from the main process.
 */

function strToIntCoords(coords)
{
     coords=coords.split(",")
     const centerX=parseFloat(coords[0])
     const centerY=parseFloat(coords[1])
     return [centerX,centerY]
}
//drawing the molecule
function drawAtom(atom,coords,ctx){
    //draws an atom given its symbol and coordinates
     centerX=strToIntCoords(coords)[0]
     centerY= strToIntCoords(coords)[1]
     ctx.save()
     ctx.fillStyle="black"
     ctx.font="40px Arial"
     ctx.textAlign="center"
     ctx.textBaseline="middle"
     ctx.fillText(atom,centerX,centerY)
     ctx.restore()
}

function drawAllAtoms(allAtoms,ctx){
    //draws all atoms of a given molecule
    for (const key in allAtoms)
    {
        drawAtom(allAtoms[key],key,ctx)
    }
}
function drawBond(type, orientation,coords,ctx)
{   
    newCoords=strToIntCoords(coords)
    centerX=newCoords[0]
    centerY=newCoords[1]
    const text=type
    ctx.save();  
    ctx.fillStyle="black"
    ctx.font="40px Arial"
    ctx.textAlign="center"
    ctx.textBaseline="middle"
    if (orientation==="|")
    {
        ctx.translate(centerX, centerY);   
        ctx.rotate(Math.PI/2)
        ctx.fillText(text,0,0)
        ctx.restore()
        return
    }
     ctx.fillText(text,centerX,centerY)
     ctx.restore()
}

function drawAllBonds(allBonds,ctx){
    for (coord in allBonds)
    {
        drawBond(allBonds[coord][0],allBonds[coord][1],coord,ctx)
    }
}


function drawAllLonePairs(allLonePairs,ctx)
{
       for (coord in allLonePairs)
       {
        drawLonePair(allLonePairs[coord][0],allLonePairs[coord][1],coord,ctx)
       }
}




