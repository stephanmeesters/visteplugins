varying vec3 Normal; 
varying vec4 EyePosition; 
uniform int NumLights; 

varying vec2 texture_coordinate;
uniform sampler2D my_color_texture;

/* Compute the contribution from a particular light source. This basically * comes straight out of the OpenGL orange book. */ 
void DirectionalLight(in int lightIndex, in vec3 normal, inout vec4 ambient, inout vec4 diffuse, inout vec4 specular) 
{ 
	/**** Compute ambient term. ****/ 
	ambient += gl_LightSource[lightIndex].ambient; 

	/**** Compute diffuse term. ****/ 
	/* normal dot light direction. Assume the light direction vector is already normalized.*/ 
	float nDotL = max(0.0, dot(normal, normalize(vec3(gl_LightSource[lightIndex].position)))); 
	diffuse += gl_LightSource[lightIndex].diffuse * nDotL; 

	/**** Compute specular term. ****/ 
	/* normal dot halfway vector */ 
	float nDotH = max(0.0, dot(normal, vec3(gl_LightSource[lightIndex].halfVector))); 

	float pf; 
	/* Power factor. */ 
	if (nDotH <= 0.0) { pf = 0.0; } 
	else { pf = pow(nDotH, gl_FrontMaterial.shininess); } 

	specular += gl_LightSource[lightIndex].specular * pf; 
} 

void AllLights(in vec3 normal, inout vec4 ambient, inout vec4 diffuse, inout vec4 specular) 
{ 
	DirectionalLight(0, normal, ambient, diffuse, specular); 
	if (NumLights > 1) 
	{ 
		DirectionalLight(1, normal, ambient, diffuse, specular);
		if (NumLights > 2) 
		{ 
			DirectionalLight(2, normal, ambient, diffuse, specular); 
			if (NumLights > 3) 
			{ 
				DirectionalLight(3, normal, ambient, diffuse, specular); 
				if (NumLights > 4) 
				{ 
					DirectionalLight(4, normal, ambient, diffuse, specular); 
				} 
			} 
		} 
	} 
} 
 
void main(void) 
{ 
	/* If lighting the back of a polygon, flip normal.*/ 
	vec3 normal = normalize(Normal); 

	/* Compute light contributions. */ 
	vec4 ambient = vec4(0.0); 
	vec4 diffuse = vec4(0.0); 
	vec4 specular = vec4(0.0); 

	AllLights(normal, ambient, diffuse, specular); 

	gl_FragColor = ( ambient*gl_FrontMaterial.ambient + diffuse*gl_FrontMaterial.diffuse + specular*gl_FrontMaterial.specular + texture2D(my_color_texture, texture_coordinate)); 
}
