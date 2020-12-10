/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://lammps.sandia.gov/, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

/* ----------------------------------------------------------------------
   Contributing authors: Axel Kohlmeyer (Temple U),
                         Ryan S. Elliott (UMN),
                         Yaser Afshar (UMN)
------------------------------------------------------------------------- */

/* ----------------------------------------------------------------------
   This program is free software; you can redistribute it and/or modify it
   under the terms of the GNU General Public License as published by the Free
   Software Foundation; either version 2 of the License, or (at your option)
   any later version.

   This program is distributed in the hope that it will be useful, but WITHOUT
   ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
   FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for
   more details.

   You should have received a copy of the GNU General Public License along with
   this program; if not, see <https://www.gnu.org/licenses>.

   Linking LAMMPS statically or dynamically with other modules is making a
   combined work based on LAMMPS. Thus, the terms and conditions of the GNU
   General Public License cover the whole combination.

   In addition, as a special exception, the copyright holders of LAMMPS give
   you permission to combine LAMMPS with free software programs or libraries
   that are released under the GNU LGPL and with code included in the standard
   release of the "kim-api" under the CDDL (or modified versions of such code,
   with unchanged license). You may copy and distribute such a system following
   the terms of the GNU GPL for LAMMPS and the licenses of the other code
   concerned, provided that you include the source code of that other code
   when and as the GNU GPL requires distribution of source code.

   Note that people who make modified versions of LAMMPS are not obligated to
   grant this special exception for their modified versions; it is their choice
   whether to do so. The GNU General Public License gives permission to release
   a modified version without this exception; this exception also makes it
   possible to release a modified version which carries forward this exception.
------------------------------------------------------------------------- */

/* ----------------------------------------------------------------------
   Designed for use with the kim-api-2.1.0 (and newer) package
------------------------------------------------------------------------- */

#include "kim_query.h"

#include "comm.h"
#include "error.h"
#include "fix_store_kim.h"
#include "info.h"
#include "input.h"
#include "modify.h"
#include "utils.h"
#include "variable.h"
#include "version.h"
#include "tokenizer.h"
#include "fmt/format.h"

#include <cstdlib>
#include <cstring>
#include <string>

#if defined(LMP_KIM_CURL)
#include <sys/types.h>
#include <curl/curl.h>
#endif

using namespace LAMMPS_NS;

#if defined(LMP_KIM_CURL)

struct WriteBuf {
  char *dataptr;
  size_t sizeleft;
};

static char *do_query(char *, char *, int, char **, int, MPI_Comm);
static size_t write_callback(void *, size_t, size_t, void *);

#endif

/* ---------------------------------------------------------------------- */

void KimQuery::command(int narg, char **arg)
{
  if (narg < 2) error->all(FLERR,"Illegal kim_query command");

  // check if we had a kim_init command by finding fix STORE/KIM
  // retrieve model name.
  char *model_name;

  const int ifix = modify->find_fix("KIM_MODEL_STORE");
  if (ifix >= 0) {
    FixStoreKIM *fix_store = (FixStoreKIM *) modify->fix[ifix];
    model_name = (char *)fix_store->getptr("model_name");
  } else error->all(FLERR,"Must use 'kim_init' before 'kim_query'");

  char *varname = arg[0];

  bool split = false;
  if (strcmp("split",arg[1]) == 0) {
    if (narg == 2) error->all(FLERR,"Illegal kim_query command.\nThe keyword "
                                    "'split' must be followed by the name of "
                                    "the query function");
    if (strcmp("list",arg[2]) == 0)
      error->all(FLERR,"Illegal kim_query command.\nThe 'list' keyword "
                       "can not be used after 'split'");
    split = true;
    arg++;
    narg--;
  }

  // The “list” is the default setting
  // the result is returned as a space-separated list of values in variable
  if (strcmp("list",arg[1]) == 0) {
    if (narg == 2) error->all(FLERR,"Illegal kim_query command.\nThe 'list' "
                                    "keyword must be followed by ('split' "
                                    "and) the name of the query function");
    arg++;
    narg--;
  }

  char *function = arg[1];
  for (int i = 2; i < narg; ++i) {
    if (strncmp("model=",arg[i],6) == 0)
      error->all(FLERR,"Illegal 'model' key in kim_query command");

    if (!strchr(arg[i], '=') || !strchr(arg[i], '[') || !strchr(arg[i], ']'))
      error->all(FLERR,fmt::format("Illegal query format.\nInput argument of "
                                   "`{}` to kim_query is wrong. The query "
                                   "format is the keyword=[value], where value "
                                   "is always an array of one or more "
                                   "comma-separated items", arg[i]));
  }

#if defined(LMP_KIM_CURL)

  char *value = do_query(function, model_name, narg-2, arg+2, comm->me, world);

  // check for valid result
  // on error the content of "value" is a '\0' byte
  // as the first element, and then the error message
  // that was returned by the web server

  if (strlen(value) == 0) {
    error->all(FLERR,fmt::format("OpenKIM query failed: {}", value+1));
  } else if (strcmp(value,"EMPTY") == 0) {
    error->all(FLERR,fmt::format("OpenKIM query returned no results"));
  }

  input->write_echo("#=== BEGIN kim-query =========================================\n");
  ValueTokenizer values(value, ",");
  if (split) {
    int counter = 1;
    while (values.has_next()) {
      auto svalue = values.next_string();
      auto setcmd = fmt::format("{}_{} string {}", varname, counter++, svalue);
      input->variable->set(setcmd);
      input->write_echo(fmt::format("variable {}\n", setcmd));
    }
  } else {
    auto svalue = values.next_string();
    std::string setcmd = fmt::format("{} string \"{}", varname, svalue);
    while (values.has_next()) {
      svalue = values.next_string();
      setcmd += fmt::format(" {}", svalue);
    }
    setcmd += "\"";
    input->variable->set(setcmd);
    input->write_echo(fmt::format("variable {}\n", setcmd));
  }
  input->write_echo("#=== END kim-query ===========================================\n\n");

  delete[] value;
#else
  error->all(FLERR,"Cannot use 'kim_query' command when KIM package "
                   "is compiled without support for libcurl");
#endif
}

#if defined(LMP_KIM_CURL)

// copy data to the user provided data structure, optionally in increments

size_t write_callback(void *data, size_t size, size_t nmemb, void *userp)
{
  struct WriteBuf *buf = (struct WriteBuf *)userp;

  // copy chunks into the buffer for as long as there is space left
  if (buf->sizeleft) {
    const size_t buffer_size = size * nmemb;
    const size_t copy_this_much =
      buf->sizeleft > buffer_size ? buffer_size : buf->sizeleft;

    memcpy(buf->dataptr, data, copy_this_much);

    buf->dataptr += copy_this_much;
    buf->sizeleft -= copy_this_much;

    return copy_this_much;
  }
  return 0; // done
}

char *do_query(char *qfunction, char * model_name, int narg, char **arg,
               int rank, MPI_Comm comm)
{
  char value[512];

  // run the web query from rank 0 only

  if (rank == 0) {
    // set up and clear receive buffer
    struct WriteBuf buf;
    buf.dataptr = value;
    buf.sizeleft = 511;
    memset(value,0,512);

    // create curl web query instance
    curl_global_init(CURL_GLOBAL_DEFAULT);
    CURL *handle = curl_easy_init();

    if (handle) {
      auto url = fmt::format("https://query.openkim.org/api/{}", qfunction);
      auto query = fmt::format("model=[\"{}\"]", model_name);
      for (int i = 0; i < narg; ++i) {
        ValueTokenizer values(arg[i], "=[]");
        std::string key = values.next_string();
        std::string val = values.next_string();
        std::string::size_type n = val.find(",");
        if (n == std::string::npos) {
          if (utils::is_integer(val) ||
              utils::is_double(val) ||
              (val.front() == '"' &&
               val.back() == '"')) {
            query += fmt::format("&{}", arg[i]);
          } else {
            query += fmt::format("&{}=[\"{}\"]", key, val);
          }
        } else {
          query += fmt::format("&{}=[", key);
          while (n != std::string::npos){
            std::string sval = val.substr(0, n);
            if (utils::is_integer(sval) ||
                utils::is_double(sval) ||
                (val.front() == '"' &&
                 val.back() == '"')) {
              query += fmt::format("{},", sval);
            } else {
              query += fmt::format("\"{}\",", sval);
            }
            val = val.substr(n + 1);
            n = val.find(",");
          }
          if (val.size()) query += fmt::format("\"{}\"]", val);
          else query[query.size() - 1]=']';
        }
      }

#if LMP_DEBUG_CURL
      curl_easy_setopt(handle, CURLOPT_VERBOSE, 1L);
#endif

#if LMP_NO_SSL_CHECK
      // Certificate Verification
      // by telling libcurl to not verify the peer.
      // Disable verifying SSL certificate and host name. Insecure.
      curl_easy_setopt(handle, CURLOPT_SSL_VERIFYPEER, 0L);
      curl_easy_setopt(handle, CURLOPT_SSL_VERIFYHOST, 0L);
#endif

      {
        char *env_c = std::getenv("CURL_CA_BUNDLE");
        if (env_c) {
          // Certificate Verification
          // by specifying your own CA cert path. Set the environment variable
          // CURL_CA_BUNDLE to the path of your choice.
          curl_easy_setopt(handle, CURLOPT_CAINFO, env_c);
        }
      }

      std::string user_agent = fmt::format("kim_query--LAMMPS/{} ({})",
                                           LAMMPS_VERSION, Info::get_os_info());

      curl_easy_setopt(handle, CURLOPT_USERAGENT, user_agent.c_str());
      curl_easy_setopt(handle, CURLOPT_URL, url.c_str());
      curl_easy_setopt(handle, CURLOPT_FOLLOWLOCATION, 1L);
      curl_easy_setopt(handle, CURLOPT_POSTFIELDS, query.c_str());
      curl_easy_setopt(handle, CURLOPT_WRITEFUNCTION, write_callback);
      curl_easy_setopt(handle, CURLOPT_WRITEDATA,&buf);

      // perform OpenKIM query and check for errors
      CURLcode res = curl_easy_perform(handle);
      if (res != CURLE_OK) {
        // on error we return an "empty" string but add error message after it
        value[0] = '\0';
        strcpy(value+1,curl_easy_strerror(res));
      }
      curl_easy_cleanup(handle);
    }
    curl_global_cleanup();
  }
  MPI_Bcast(value, 512, MPI_CHAR, 0, comm);

  // we must make a proper copy of the query, as the stack allocation
  // for "value" will go out of scope. a valid query has a '[' as
  // the first character. skip over it (and the last character ']', too)
  // an error message starts with a '\0' character. copy that and
  // the remaining string, as that is the error message.

  char *retval;
  // a valid query has a '[' as the first character
  if (value[0] == '[') {
    int len = strlen(value) - 1;
    if (value[len] == ']') {
      value[len] = '\0';
      retval = new char[len];
      if (strcmp(value+1, "") == 0) strcpy(retval,"EMPTY");
      else strcpy(retval,value+1);
    } else {
      retval = new char[len+2];
      retval[0] = '\0';
      strcpy(retval+1,value);
    }
  // an error message starts with a '\0' character
  } else if (value[0] == '\0') {
    int len = strlen(value+1)+2;
    retval = new char[len];
    retval[0] = '\0';
    strcpy(retval+1,value+1);
  // unknown response type. we should not get here.
  } else {
    // we return an "empty" string but add error message after it
    int len = strlen(value)+2;
    retval = new char[len];
    retval[0] = '\0';
    strcpy(retval+1,value);
  }
  return retval;
}
#endif
