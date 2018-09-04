#include "Game.hpp"

#include "gl_errors.hpp" //helper for dumpping OpenGL error messages
#include "read_chunk.hpp" //helper for reading a vector of structures from a file
#include "data_path.hpp" //helper to get paths relative to executable

#include <glm/gtc/type_ptr.hpp>

#include <iostream>
#include <fstream>
#include <map>
#include <cstddef>
#include <random>

//helper defined later; throws if shader compilation fails:
static GLuint compile_shader(GLenum type, std::string const &source);

Game::Game() {
	{ //create an opengl program to perform sun/sky (well, directional+hemispherical) lighting:
		GLuint vertex_shader = compile_shader(GL_VERTEX_SHADER,
			"#version 330\n"
			"uniform mat4 object_to_clip;\n"
			"uniform mat4x3 object_to_light;\n"
			"uniform mat3 normal_to_light;\n"
			"layout(location=0) in vec4 Position;\n" //note: layout keyword used to make sure that the location-0 attribute is always bound to something
			"in vec3 Normal;\n"
			"in vec4 Color;\n"
			"out vec3 position;\n"
			"out vec3 normal;\n"
			"out vec4 color;\n"
			"void main() {\n"
			"	gl_Position = object_to_clip * Position;\n"
			"	position = object_to_light * Position;\n"
			"	normal = normal_to_light * Normal;\n"
			"	color = Color;\n"
			"}\n"
		);

		GLuint fragment_shader = compile_shader(GL_FRAGMENT_SHADER,
			"#version 330\n"
			"uniform vec3 sun_direction;\n"
			"uniform vec3 sun_color;\n"
			"uniform vec3 sky_direction;\n"
			"uniform vec3 sky_color;\n"
			"in vec3 position;\n"
			"in vec3 normal;\n"
			"in vec4 color;\n"
			"out vec4 fragColor;\n"
			"void main() {\n"
			"	vec3 total_light = vec3(0.0, 0.0, 0.0);\n"
			"	vec3 n = normalize(normal);\n"
			"	{ //sky (hemisphere) light:\n"
			"		vec3 l = sky_direction;\n"
			"		float nl = 0.5 + 0.5 * dot(n,l);\n"
			"		total_light += nl * sky_color;\n"
			"	}\n"
			"	{ //sun (directional) light:\n"
			"		vec3 l = sun_direction;\n"
			"		float nl = max(0.0, dot(n,l));\n"
			"		total_light += nl * sun_color;\n"
			"	}\n"
			"	fragColor = vec4(color.rgb * total_light, color.a);\n"
			"}\n"
		);

		simple_shading.program = glCreateProgram();
		glAttachShader(simple_shading.program, vertex_shader);
		glAttachShader(simple_shading.program, fragment_shader);
		//shaders are reference counted so this makes sure they are freed after program is deleted:
		glDeleteShader(vertex_shader);
		glDeleteShader(fragment_shader);

		//link the shader program and throw errors if linking fails:
		glLinkProgram(simple_shading.program);
		GLint link_status = GL_FALSE;
		glGetProgramiv(simple_shading.program, GL_LINK_STATUS, &link_status);
		if (link_status != GL_TRUE) {
			std::cerr << "Failed to link shader program." << std::endl;
			GLint info_log_length = 0;
			glGetProgramiv(simple_shading.program, GL_INFO_LOG_LENGTH, &info_log_length);
			std::vector< GLchar > info_log(info_log_length, 0);
			GLsizei length = 0;
			glGetProgramInfoLog(simple_shading.program, GLsizei(info_log.size()), &length, &info_log[0]);
			std::cerr << "Info log: " << std::string(info_log.begin(), info_log.begin() + length);
			throw std::runtime_error("failed to link program");
		}
	}

	{ //read back uniform and attribute locations from the shader program:
		simple_shading.object_to_clip_mat4 = glGetUniformLocation(simple_shading.program, "object_to_clip");
		simple_shading.object_to_light_mat4x3 = glGetUniformLocation(simple_shading.program, "object_to_light");
		simple_shading.normal_to_light_mat3 = glGetUniformLocation(simple_shading.program, "normal_to_light");

		simple_shading.sun_direction_vec3 = glGetUniformLocation(simple_shading.program, "sun_direction");
		simple_shading.sun_color_vec3 = glGetUniformLocation(simple_shading.program, "sun_color");
		simple_shading.sky_direction_vec3 = glGetUniformLocation(simple_shading.program, "sky_direction");
		simple_shading.sky_color_vec3 = glGetUniformLocation(simple_shading.program, "sky_color");

		simple_shading.Position_vec4 = glGetAttribLocation(simple_shading.program, "Position");
		simple_shading.Normal_vec3 = glGetAttribLocation(simple_shading.program, "Normal");
		simple_shading.Color_vec4 = glGetAttribLocation(simple_shading.program, "Color");
	}

	struct Vertex {
		glm::vec3 Position;
		glm::vec3 Normal;
		glm::u8vec4 Color;
	};
	static_assert(sizeof(Vertex) == 28, "Vertex should be packed.");

	{ //load mesh data from a binary blob:
		std::ifstream blob(data_path("meshes.blob"), std::ios::binary);
		//The blob will be made up of three chunks:
		// the first chunk will be vertex data (interleaved position/normal/color)
		// the second chunk will be characters
		// the third chunk will be an index, mapping a name (range of characters) to a mesh (range of vertex data)

		//read vertex data:
		std::vector< Vertex > vertices;
		read_chunk(blob, "dat0", &vertices);

		//read character data (for names):
		std::vector< char > names;
		read_chunk(blob, "str0", &names);

		//read index:
		struct IndexEntry {
			uint32_t name_begin;
			uint32_t name_end;
			uint32_t vertex_begin;
			uint32_t vertex_end;
		};
		static_assert(sizeof(IndexEntry) == 16, "IndexEntry should be packed.");

		std::vector< IndexEntry > index_entries;
		read_chunk(blob, "idx0", &index_entries);

		if (blob.peek() != EOF) {
			std::cerr << "WARNING: trailing data in meshes file." << std::endl;
		}

		//upload vertex data to the graphics card:
		glGenBuffers(1, &meshes_vbo);
		glBindBuffer(GL_ARRAY_BUFFER, meshes_vbo);
		glBufferData(GL_ARRAY_BUFFER, sizeof(Vertex) * vertices.size(), vertices.data(), GL_STATIC_DRAW);
		glBindBuffer(GL_ARRAY_BUFFER, 0);

		//create map to store index entries:
		std::map< std::string, Mesh > index;
		for (IndexEntry const &e : index_entries) {
			if (e.name_begin > e.name_end || e.name_end > names.size()) {
				throw std::runtime_error("invalid name indices in index.");
			}
			if (e.vertex_begin > e.vertex_end || e.vertex_end > vertices.size()) {
				throw std::runtime_error("invalid vertex indices in index.");
			}
			Mesh mesh;
			mesh.first = e.vertex_begin;
			mesh.count = e.vertex_end - e.vertex_begin;
			auto ret = index.insert(std::make_pair(
				std::string(names.begin() + e.name_begin, names.begin() + e.name_end),
				mesh));
			if (!ret.second) {
				throw std::runtime_error("duplicate name in index.");
			}
		}

		//look up into index map to extract meshes:
		auto lookup = [&index](std::string const &name) -> Mesh {
			auto f = index.find(name);
			if (f == index.end()) {
				throw std::runtime_error("Mesh named '" + name + "' does not appear in index.");
			}
			return f->second;
		};

		floor_mesh = lookup("Floor");
		plus_mesh = lookup("+");
		L_mesh = lookup("L");
		I_mesh = lookup("I");
		i_mesh = lookup("i");
		T_mesh = lookup("T");
		player_mesh = lookup("Player");
		goal_mesh = lookup("Goal");
	}

	{ //create vertex array object to hold the map from the mesh vertex buffer to shader program attributes:
		glGenVertexArrays(1, &meshes_for_simple_shading_vao);
		glBindVertexArray(meshes_for_simple_shading_vao);
		glBindBuffer(GL_ARRAY_BUFFER, meshes_vbo);
		//note that I'm specifying a 3-vector for a 4-vector attribute here, and this is okay to do:
		glVertexAttribPointer(simple_shading.Position_vec4, 3, GL_FLOAT, GL_FALSE, sizeof(Vertex), (GLbyte *)0 + offsetof(Vertex, Position));
		glEnableVertexAttribArray(simple_shading.Position_vec4);
		if (simple_shading.Normal_vec3 != -1U) {
			glVertexAttribPointer(simple_shading.Normal_vec3, 3, GL_FLOAT, GL_FALSE, sizeof(Vertex), (GLbyte *)0 + offsetof(Vertex, Normal));
			glEnableVertexAttribArray(simple_shading.Normal_vec3);
		}
		if (simple_shading.Color_vec4 != -1U) {
			glVertexAttribPointer(simple_shading.Color_vec4, 4, GL_UNSIGNED_BYTE, GL_TRUE, sizeof(Vertex), (GLbyte *)0 + offsetof(Vertex, Color));
			glEnableVertexAttribArray(simple_shading.Color_vec4);
		}
		glBindBuffer(GL_ARRAY_BUFFER, 0);
	}

	GL_ERRORS();

	//----------------

	// by default, slight shear left
	shear = glm::angleAxis(0.3f, glm::vec3(-0.25f, 0.25f, 0.0f));
	board_rotation = glm::normalize(shear * glm::quat(0.0f, 0.0f, 1.0f, 0.0f));

	generate_board();
}

Game::~Game() {
	glDeleteVertexArrays(1, &meshes_for_simple_shading_vao);
	meshes_for_simple_shading_vao = -1U;

	glDeleteBuffers(1, &meshes_vbo);
	meshes_vbo = -1U;

	glDeleteProgram(simple_shading.program);
	simple_shading.program = -1U;

	GL_ERRORS();
}

void Game::set_mesh(glm::vec2 pos, int connected_dirs_bitvec) {
	bool up = get_dir_bitvec(connected_dirs_bitvec, Game::UP);
	bool left = get_dir_bitvec(connected_dirs_bitvec, Game::LEFT);
	bool down = get_dir_bitvec(connected_dirs_bitvec, Game::DOWN);
	bool right = get_dir_bitvec(connected_dirs_bitvec, Game::RIGHT);

	// update connectivity
	tile_connections[pos.y*board_size.x+pos.x] = connected_dirs_bitvec;

	// update meshes
	if (up && down && !left && !right) { // I
		path_meshes[pos.y*board_size.x+pos.x] = &I_mesh;
	}
	else if (!up && !down && left && right) { // I
		path_meshes[pos.y*board_size.x+pos.x] = &I_mesh;
		glm::quat &r = tile_rotations[pos.y*board_size.x+pos.x];
		r = glm::normalize(r * glm::angleAxis(0.5f*glm::pi<float>(), glm::vec3(0.0f, 0.0f, 1.0f)));
	}
	else if (((int)up + (int)left + (int)down + (int)right) == 1) { // i
		path_meshes[pos.y*board_size.x+pos.x] = &i_mesh;
		glm::quat &r = tile_rotations[pos.y*board_size.x+pos.x];
		int rot = (int)down * (int)Game::UP + (int)left * (int)Game::LEFT + (int)up * (int)Game::DOWN + (int)right * Game::RIGHT;
		r = glm::normalize(r * glm::angleAxis(rot*0.5f*glm::pi<float>(), glm::vec3(0.0f, 0.0f, 1.0f)));
	}
	else if (!up && down && left && right) { // T
		path_meshes[pos.y*board_size.x+pos.x] = &T_mesh;
	}
	else if (up && down && left && !right) { // T rotated 90
		path_meshes[pos.y*board_size.x+pos.x] = &T_mesh;
		glm::quat &r = tile_rotations[pos.y*board_size.x+pos.x];
		r = glm::normalize(r * glm::angleAxis(0.5f*glm::pi<float>(), glm::vec3(0.0f, 0.0f, 1.0f)));
	}
	else if (up && !down && left && right) { // T rotated 180
		path_meshes[pos.y*board_size.x+pos.x] = &T_mesh;
		glm::quat &r = tile_rotations[pos.y*board_size.x+pos.x];
		r = glm::normalize(r * glm::angleAxis(glm::pi<float>(), glm::vec3(0.0f, 0.0f, 1.0f)));
	}
	else if (up && down && !left && right) { // T rotated 270
		path_meshes[pos.y*board_size.x+pos.x] = &T_mesh;
		glm::quat &r = tile_rotations[pos.y*board_size.x+pos.x];
		r = glm::normalize(r * glm::angleAxis(1.5f*glm::pi<float>(), glm::vec3(0.0f, 0.0f, 1.0f)));
	}
	else if (!up && down && !left && right) { // r
		path_meshes[pos.y*board_size.x+pos.x] = &L_mesh;
		glm::quat &r = tile_rotations[pos.y*board_size.x+pos.x];
		r = glm::normalize(r * glm::angleAxis(1.5f*glm::pi<float>(), glm::vec3(0.0f, 0.0f, 1.0f)));
	}
	else if (!up && down && left && !right) { // r rotated 90
		path_meshes[pos.y*board_size.x+pos.x] = &L_mesh;
	}
	else if (up && !down && left && !right) { // r rotated 180
		path_meshes[pos.y*board_size.x+pos.x] = &L_mesh;
		glm::quat &r = tile_rotations[pos.y*board_size.x+pos.x];
		r = glm::normalize(r * glm::angleAxis(0.5f*glm::pi<float>(), glm::vec3(0.0f, 0.0f, 1.0f)));
	}
	else if (up && !down && !left && right) { // r rotated 270
		path_meshes[pos.y*board_size.x+pos.x] = &L_mesh;
		glm::quat &r = tile_rotations[pos.y*board_size.x+pos.x];
		r = glm::normalize(r * glm::angleAxis(glm::pi<float>(), glm::vec3(0.0f, 0.0f, 1.0f)));
	}
	else { // +
		path_meshes[pos.y*board_size.x+pos.x] = &plus_mesh;
	}
}

bool Game::valid_pos(glm::vec2 pos) {
	return (0 <= pos.x && pos.x < board_size.x) && (0 <= pos.y && pos.y < board_size.y);
}

bool Game::equal_pos(glm::vec2 first, glm::vec2 second) {
	return (int)first.x == (int)second.x && (int)first.y == (int)second.y;
}

void Game::set_dir_bitvec(int &dir_bitvec, DIR dir) {
	dir_bitvec |= (1 << dir);
}
bool Game::get_dir_bitvec(int dir_bitvec, DIR dir) {
	return dir_bitvec & (1 << dir);
}

bool Game::generate_board_helper(glm::vec2 pos, glm::vec2 from_dir, std::unordered_set<int> &visited) {
	// mark as visited
	visited.insert(pos.y*board_size.x+pos.x);

	// std::cout << "visiting " << pos.x << " " << pos.y << ", size: " << visited.size() << std::endl;

	// keep track of whether dfs completed or not
	bool done = false;

	// keep track of which directions are connected
	int connected_dirs_bitvec = 0;

	// the moment we visit the last cell, set the cell as goal, and backtrack to
	// the starting cell. every cell traversed during this final "backtrack" will
	// be included in the resulting path.
	if (visited.size() >= (board_size.x * board_size.y)) {
		goal = pos;
		done = true;
	}

	// generate permuted list of directions
	std::vector<DIR> shuffled_dirs = {UP,LEFT,DOWN,RIGHT};
	std::random_shuffle(shuffled_dirs.begin(), shuffled_dirs.end());

	for (auto it = shuffled_dirs.begin(); it != shuffled_dirs.end(); it++) {
		glm::vec2 dir = direction_vecs[*it];
		bool through = false;

		// std::cout << dir.x << "," << dir.y << " " << from_dir.x << "," << from_dir.y << std::endl;

		if (equal_pos(dir, -from_dir)) {
			through = true;
		} else {
			glm::vec2 new_pos = pos + dir;
			if (visited.find(new_pos.y*board_size.x+new_pos.x) == visited.end() && valid_pos(new_pos)) {
				if (generate_board_helper(new_pos, dir, visited)) {
					through = true;
					done = true;
				}
			}
		}

		if (through) {
			if (dir.y > 0) { set_dir_bitvec(connected_dirs_bitvec, Game::UP); }
			if (dir.x < 0) { set_dir_bitvec(connected_dirs_bitvec, Game::LEFT); }
			if (dir.y < 0) { set_dir_bitvec(connected_dirs_bitvec, Game::DOWN); }
			if (dir.x > 0) { set_dir_bitvec(connected_dirs_bitvec, Game::RIGHT); }
		}
	}

	if (done) {
		// std::cout << pos.x << " " << pos.y << std::endl;
		set_mesh(pos, connected_dirs_bitvec);
	}
	return done;
}

void Game::clear_board() {
	path_meshes.clear();
	tile_rotations.clear();
	tile_connections.clear();

	// re-reserve to maintain capacity
	path_meshes.reserve(board_size.x * board_size.y);
	tile_rotations.reserve(board_size.x * board_size.y);
	tile_connections.reserve(board_size.x * board_size.y);

	for (uint32_t i = 0; i < board_size.x * board_size.y; ++i) {
		path_meshes.emplace_back(&floor_mesh);
		tile_rotations.emplace_back(board_rotation);
		tile_connections.emplace_back(0);
	}
}

void Game::generate_board() {
	std::mt19937 mt(time(NULL));

	// pick a random position on the bottom or top edge of board
	float x = mt()%board_size.x;
	float y = mt()%2 * (board_size.y-1);

	start = glm::vec2(x, y);
	player = start;

	clear_board();

	// scramble keys
	std::random_shuffle(shuffled_keys.begin(), shuffled_keys.end());

	std::unordered_set<int> visited;
	generate_board_helper(player, glm::vec2(0,0), visited);
}

Game::DIR Game::flip_dir(Game::DIR dir) {
	return (Game::DIR)((dir + 2) % 4);
}

bool Game::handle_event(SDL_Event const &evt, glm::uvec2 window_size) {
	//ignore any keys that are the result of automatic key repeat:
	if ((evt.type == SDL_KEYDOWN && evt.key.repeat) || reset_after_updates > 0) {
		return false;
	}
	//move cursor on L/R/U/D press:
	if (evt.type == SDL_KEYDOWN && evt.key.repeat == 0) {
		DIR dir = Game::STAY;
		if (evt.key.keysym.scancode == SDL_SCANCODE_UP) {
			dir = shuffled_keys[0];
		} else if (evt.key.keysym.scancode == SDL_SCANCODE_LEFT) {
			dir = shuffled_keys[1];
		} else if (evt.key.keysym.scancode == SDL_SCANCODE_DOWN) {
			dir = shuffled_keys[2];
		} else if (evt.key.keysym.scancode == SDL_SCANCODE_RIGHT) {
			dir = shuffled_keys[3];
		}

		if (dir == Game::STAY || !valid_pos(player + direction_vecs[dir])) {
			return false;
		}

		if (!get_dir_bitvec(tile_connections[player.y*board_size.x+player.x], dir)) {
			reset_after_updates = 20;
			player += (direction_vecs[dir] / 2.0f);
		} else {
			player += (direction_vecs[dir]);
		}
		return true;
	}

	return false;
}

void Game::update(float elapsed) {
	if (reset_after_updates > 0) {
		reset_after_updates--;

		if (reset_after_updates == 0) {
			// reset player to start
			player = start;
		}
	}

	if (player == goal) {
		score++;
		// generate new level
		generate_board();
	}
}

void Game::draw(glm::uvec2 drawable_size) {
	//Set up a transformation matrix to fit the board in the window:
	glm::mat4 world_to_clip;
	{
		float aspect = float(drawable_size.x) / float(drawable_size.y);

		//want scale such that board * scale fits in [-aspect,aspect]x[-1.0,1.0] screen box:
		float scale = glm::min(
			2.0f * aspect / float(board_size.x),
			2.0f / float(board_size.y)
		);

		//center of board will be placed at center of screen:
		glm::vec2 center = 0.5f * glm::vec2(board_size);

		//NOTE: glm matrices are specified in column-major order
		world_to_clip = glm::mat4(
			scale / aspect, 0.0f, 0.0f, 0.0f,
			0.0f, scale, 0.0f, 0.0f,
			0.0f, 0.0f,-1.0f, 0.0f,
			-(scale / aspect) * center.x, -scale * center.y, 0.0f, 1.0f
		);
	}

	//set up graphics pipeline to use data from the meshes and the simple shading program:
	glBindVertexArray(meshes_for_simple_shading_vao);
	glUseProgram(simple_shading.program);

	glUniform3fv(simple_shading.sun_color_vec3, 1, glm::value_ptr(glm::vec3(0.81f, 0.81f, 0.76f)));
	glUniform3fv(simple_shading.sun_direction_vec3, 1, glm::value_ptr(glm::normalize(glm::vec3(-0.2f, 0.2f, 1.0f))));
	glUniform3fv(simple_shading.sky_color_vec3, 1, glm::value_ptr(glm::vec3(0.2f, 0.2f, 0.3f)));
	glUniform3fv(simple_shading.sky_direction_vec3, 1, glm::value_ptr(glm::vec3(0.0f, 1.0f, 0.0f)));

	//helper function to draw a given mesh with a given transformation:
	auto draw_mesh = [&](Mesh const &mesh, glm::mat4 const &object_to_world) {
		//set up the matrix uniforms:
		if (simple_shading.object_to_clip_mat4 != -1U) {
			glm::mat4 object_to_clip = world_to_clip * object_to_world;
			glUniformMatrix4fv(simple_shading.object_to_clip_mat4, 1, GL_FALSE, glm::value_ptr(object_to_clip));
		}
		if (simple_shading.object_to_light_mat4x3 != -1U) {
			glUniformMatrix4x3fv(simple_shading.object_to_light_mat4x3, 1, GL_FALSE, glm::value_ptr(object_to_world));
		}
		if (simple_shading.normal_to_light_mat3 != -1U) {
			//NOTE: if there isn't any non-uniform scaling in the object_to_world matrix, then the inverse transpose is the matrix itself, and computing it wastes some CPU time:
			glm::mat3 normal_to_world = glm::inverse(glm::transpose(glm::mat3(object_to_world)));
			glUniformMatrix3fv(simple_shading.normal_to_light_mat3, 1, GL_FALSE, glm::value_ptr(normal_to_world));
		}

		//draw the mesh:
		glDrawArrays(GL_TRIANGLES, mesh.first, mesh.count);
	};

	for (uint32_t y = 0; y < board_size.y; ++y) {
		for (uint32_t x = 0; x < board_size.x; ++x) {
			draw_mesh(floor_mesh,
				glm::mat4(
					1.0f, 0.0f, 0.0f, 0.0f,
					0.0f, 1.0f, 0.0f, 0.0f,
					0.0f, 0.0f, 1.0f, 0.0f,
					x+0.5f, y+0.5f,-0.5f, 1.0f
				)
			);
			if (path_meshes[y*board_size.x+x] != NULL) {
				draw_mesh(*path_meshes[y*board_size.x+x],
					glm::mat4(
						0.5f, 0.0f, 0.0f, 0.0f,
						0.0f, 0.5f, 0.0f, 0.0f,
						0.0f, 0.0f, 0.5f, 0.0f,
						x+0.5f, y+0.5f, 0.0f, 1.0f
					)
					* glm::mat4_cast(glm::normalize(tile_rotations[y*board_size.x+x]))
				);
			}
		}
	}
	draw_mesh(player_mesh,
		glm::mat4(
			0.5f, 0.0f, 0.0f, 0.0f,
			0.0f, 0.5f, 0.0f, 0.0f,
			0.0f, 0.0f, 0.5f, 0.0f,
			player.x+0.5f, player.y+0.5f, 0.0f, 1.0f
		)
		* glm::mat4_cast(shear)
	);
	draw_mesh(goal_mesh,
		glm::mat4(
			0.5f, 0.0f, 0.0f, 0.0f,
			0.0f, 0.5f, 0.0f, 0.0f,
			0.0f, 0.0f, 0.5f, 0.0f,
			goal.x+0.5f, goal.y+0.5f, 0.0f, 1.0f
		)
	);
	// display the score as stars to the side of the board; if score gets too
	// high start shrinking the stars
	// 9 is the max number that fits in its default size
	float size = std::min((9.0f / (float)score), 1.0f) * 0.5f;
	for (int i = 0; i < score; i++) {
		draw_mesh(goal_mesh,
			glm::mat4(
				size, 0.0f, 0.0f, 0.0f,
				0.0f, size, 0.0f, 0.0f,
				0.0f, 0.0f, size, 0.0f,
				-1.0f, (i+1)*size, 0.0f, 1.0f
			)
		);
	}

	glUseProgram(0);

	GL_ERRORS();
}



//create and return an OpenGL vertex shader from source:
static GLuint compile_shader(GLenum type, std::string const &source) {
	GLuint shader = glCreateShader(type);
	GLchar const *str = source.c_str();
	GLint length = GLint(source.size());
	glShaderSource(shader, 1, &str, &length);
	glCompileShader(shader);
	GLint compile_status = GL_FALSE;
	glGetShaderiv(shader, GL_COMPILE_STATUS, &compile_status);
	if (compile_status != GL_TRUE) {
		std::cerr << "Failed to compile shader." << std::endl;
		GLint info_log_length = 0;
		glGetShaderiv(shader, GL_INFO_LOG_LENGTH, &info_log_length);
		std::vector< GLchar > info_log(info_log_length, 0);
		GLsizei length = 0;
		glGetShaderInfoLog(shader, GLsizei(info_log.size()), &length, &info_log[0]);
		std::cerr << "Info log: " << std::string(info_log.begin(), info_log.begin() + length);
		glDeleteShader(shader);
		throw std::runtime_error("Failed to compile shader.");
	}
	return shader;
}
