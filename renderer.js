function strToIntCoords(coords)
{
     coords=coords.split(",")
     const centerX=parseFloat(coords[0])
     const centerY=parseFloat(coords[1])
     return [centerX,centerY]
}
function drawCircle(ctx, x, y, radius, color) {
    ctx.beginPath();               // start a new path
    ctx.arc(x, y, radius, 0, Math.PI * 2); // full circle
    ctx.fillStyle = color;         // fill color
    ctx.fill();                    // fill the circle
}
function getDistance(x1, y1, x2, y2) {
    const dx = x2 - x1;
    const dy = y2 - y1;
    return Math.sqrt(dx * dx + dy * dy);
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
    ctx.font="40px 'JetBrains Mono'"
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
function drawInfo(info, ctx) {
    ctx.save();

    const lineHeight = 22;
    const padding = 12;
    const startX = ctx.canvas.width - 10; // right margin
    const startY = 10;

    ctx.font = "16px 'JetBrains Mono', monospace";
    ctx.textAlign = "right";
    ctx.textBaseline = "top";

    // Prepare text lines
    const lines = Object.keys(info).map(key => `${key}: ${info[key]}`);
    const textWidth = Math.max(...lines.map(text => ctx.measureText(text).width));
    const rectWidth = textWidth + padding * 2;
    const rectHeight = lines.length * lineHeight + padding * 2;

    // Draw semi-transparent rounded rectangle WITHOUT border
    const radius = 10;
    ctx.fillStyle = "rgba(255, 255, 255, 0.9)";
    ctx.shadowColor = "transparent";
    ctx.shadowBlur = 8;
    ctx.shadowOffsetX = 2;
    ctx.shadowOffsetY = 2;

    ctx.beginPath();
    ctx.moveTo(startX - rectWidth + radius, startY);
    ctx.lineTo(startX - radius, startY);
    ctx.quadraticCurveTo(startX, startY, startX, startY + radius);
    ctx.lineTo(startX, startY + rectHeight - radius);
    ctx.quadraticCurveTo(startX, startY + rectHeight, startX - radius, startY + rectHeight);
    ctx.lineTo(startX - rectWidth + radius, startY + rectHeight);
    ctx.quadraticCurveTo(startX - rectWidth, startY + rectHeight, startX - rectWidth, startY + rectHeight - radius);
    ctx.lineTo(startX - rectWidth, startY + radius);
    ctx.quadraticCurveTo(startX - rectWidth, startY, startX - rectWidth + radius, startY);
    ctx.closePath();

    ctx.fill();  // <-- fill only, no stroke

    // Draw text
    ctx.fillStyle = "#333";
    let yPos = startY + padding;
    for (const text of lines) {
        ctx.fillText(text, startX - padding, yPos);
        yPos += lineHeight;
    }

    ctx.restore();

    return { x: startX - rectWidth, y: startY, width: rectWidth, height: rectHeight };
}


// Clear the info panel
function clearInfoPanel(ctx, rect) {
    ctx.clearRect(rect.x, rect.y, rect.width, rect.height);
}

function drawHelpButton(ctx) {
    const x = ctx.canvas.width - 30; // center x
    const y = 30;                     // center y
    const radius = 20;

    // Draw circle
    ctx.save();
    ctx.fillStyle = "blue";
    ctx.shadowColor = "rgba(0,0,0,0.3)";
    ctx.shadowBlur = 5;
    ctx.shadowOffsetX = 2;
    ctx.shadowOffsetY = 2;
    ctx.beginPath();
    ctx.arc(x, y, radius, 0, Math.PI * 2);
    ctx.fill();

    // Draw question mark
    ctx.shadowColor = "transparent"; // no shadow for text
    ctx.fillStyle = "white";
    ctx.font = "20px 'JetBrains Mono', monospace";
    ctx.textAlign = "center";
    ctx.textBaseline = "middle";
    ctx.fillText("?", x, y);

    ctx.restore();

    return { x: x, y: y, radius: radius }; // return for click detection
}


