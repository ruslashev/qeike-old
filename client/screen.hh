#pragma once

#define GLEW_STATIC
#include <GL/glew.h>
#include <SDL2/SDL.h>
#include <string>

class screen {
  SDL_Window *_window;
  SDL_GLContext _gl_context;
  std::string _title;
  int _pre_lock_mouse_x, _pre_lock_mouse_y;
public:
  int window_width, window_height;
  bool running;

  screen(const std::string &n_title, int n_window_width, int n_window_height);
  ~screen();
  void mainloop(void (*load_cb)(screen*)
      , void (*key_event_cb)(char, bool)
      , void (*mouse_motion_event_cb)(float, float, int, int)
      , void (*mouse_button_event_cb)(int, bool, int, int)
      , void (*update_cb)(double, double, screen*)
      , void (*draw_cb)(double)
      , void (*cleanup_cb)(void));
  void lock_mouse();
  void unlock_mouse();
  double get_time_in_seconds();
};

