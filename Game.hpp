#pragma once

#include "GL.hpp"

#include <SDL.h>
#include <glm/glm.hpp>
#include <glm/gtc/quaternion.hpp>

#include <vector>
#include <unordered_set>
#include <random>

// The 'Game' struct holds all of the game-relevant state,
// and is called by the main loop.

struct Game {
	//Game creates OpenGL resources (i.e. vertex buffer objects) in its
	//constructor and frees them in its destructor.
	Game();
	~Game();

	//handle_event is called when new mouse or keyboard events are received:
	// (note that this might be many times per frame or never)
	//The function should return 'true' if it handled the event.
	bool handle_event(SDL_Event const &evt, glm::uvec2 window_size);

	//update is called at the start of a new frame, after events are handled:
	void update(float elapsed);

	//draw is called after update:
	void draw(glm::uvec2 drawable_size);

	enum DIR { STAY=-1, UP, LEFT, DOWN, RIGHT };

	void set_mesh(glm::vec2 pos, int connected_dirs_bitvec);
	bool valid_pos(glm::vec2 pos);
	bool equal_pos(glm::vec2 first, glm::vec2 second);
	void set_dir_bitvec(int &dir_bitvec, DIR dir);
	bool get_dir_bitvec(int dir_bitvec, DIR dir);
	DIR flip_dir(DIR dir);
	bool generate_board_helper(glm::vec2 pos, glm::vec2 from_dir, std::unordered_set<int> &visited);
	void generate_board();
	void clear_board();

	//------- opengl resources -------

	//shader program that draws lit objects with vertex colors:
	struct {
		GLuint program = -1U; //program object

		//uniform locations:
		GLuint object_to_clip_mat4 = -1U;
		GLuint object_to_light_mat4x3 = -1U;
		GLuint normal_to_light_mat3 = -1U;
		GLuint sun_direction_vec3 = -1U;
		GLuint sun_color_vec3 = -1U;
		GLuint sky_direction_vec3 = -1U;
		GLuint sky_color_vec3 = -1U;

		//attribute locations:
		GLuint Position_vec4 = -1U;
		GLuint Normal_vec3 = -1U;
		GLuint Color_vec4 = -1U;
	} simple_shading;

	//mesh data, stored in a vertex buffer:
	GLuint meshes_vbo = -1U; //vertex buffer holding mesh data

	//The location of each mesh in the meshes vertex buffer:
	struct Mesh {
		GLint first = 0;
		GLsizei count = 0;
	};

	Mesh floor_mesh;
	Mesh plus_mesh;
	Mesh L_mesh;
	Mesh I_mesh;
	Mesh i_mesh;
	Mesh T_mesh;
	Mesh player_mesh;
	Mesh goal_mesh;

	GLuint meshes_for_simple_shading_vao = -1U; //vertex array object that describes how to connect the meshes_vbo to the simple_shading_program

	//------- game state -------

	glm::uvec2 board_size = glm::uvec2(5,5);
	glm::quat shear;
	glm::quat board_rotation;

	// up, left, down, right
	std::vector<glm::vec2> direction_vecs = {glm::vec2(0,1), glm::vec2(-1,0), glm::vec2(0,-1), glm::vec2(1,0)};

	std::vector< Mesh * > path_meshes;
	std::vector< glm::quat > tile_rotations;
	std::vector< int > tile_connections;
	std::vector< DIR > shuffled_keys = {UP,LEFT,DOWN,RIGHT};

	glm::vec2 start = glm::vec2(0,0);
	glm::vec2 player = glm::vec2(0,0);
	glm::vec2 goal = glm::vec2(0,0);

	int reset_after_updates = 0;
	int score = 0;

};
