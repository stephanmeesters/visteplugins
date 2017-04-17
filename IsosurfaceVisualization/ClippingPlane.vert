uniform float ClipX;
uniform float ClipY;
uniform float ClipZ;

varying vec3 Normal; 
varying vec4 EyePosition; 

varying vec2 texture_coordinate;

void main(void) 
{ 
	if(gl_Vertex[0] > ClipX)
		gl_ClipDistance[0] = 1;
	else
		gl_ClipDistance[0] = -1;

	if(gl_Vertex[1] > ClipY)
		gl_ClipDistance[1] = 1;
	else
		gl_ClipDistance[1] = -1;

	if(gl_Vertex[2] > ClipZ)
		gl_ClipDistance[2] = 1;
	else
		gl_ClipDistance[2] = -1;
		
	gl_Position = ftransform(); /* Compute the position in eye coordinates. */ 
	EyePosition = gl_ModelViewMatrix*gl_Vertex; /* Transform the normal. */ 
	Normal = normalize(gl_NormalMatrix*gl_Normal); /* Pass lighting and coloring parameters. */ 
	gl_FrontColor = gl_Color; 
	gl_BackColor = gl_Color; 
	gl_TexCoord[0] = gl_MultiTexCoord0; 

	texture_coordinate = vec2(gl_MultiTexCoord0);
}
