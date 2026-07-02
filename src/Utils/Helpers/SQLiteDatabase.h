#ifndef HELPERS_SQLITEDATABASE_H
#define HELPERS_SQLITEDATABASE_H

#include <sqlite3.h>
#include <spdlog/spdlog.h>

#include <ctime>
#include <string>

// NOTE: Consider using other SQLite wrapper library for better binding and
// error handling
// SQLiteDatabase class manages the lifecycle and operations of an SQLite
// database connection.
class SQLiteDatabase {
  //disallow copy and move
  SQLiteDatabase(const SQLiteDatabase&) = delete;
  void operator=(const SQLiteDatabase&) = delete;
private:
  sqlite3* db_ = nullptr;  // Pointer to the SQLite database

  // Binds an integer value to the first placeholder in the prepared SQL
  // statement.
  // stmt: Pointer to the prepared SQL statement.
  // value: Integer value to bind.
  static void bind_single_value(sqlite3_stmt* stmt, int index, int value) {
    sqlite3_bind_int(stmt, index, value);
  }

  // Binds a 16-bit unsigned integer value
  static void bind_single_value(sqlite3_stmt* stmt, int index, std::uint16_t value) {
    sqlite3_bind_int(stmt, index, static_cast<int>(value));
  }

  // Binds a 16-bit signed integer value
  static void bind_single_value(sqlite3_stmt* stmt, int index, std::int16_t value) {
    sqlite3_bind_int(stmt, index, static_cast<int>(value));
  }

  // Binds a 32-bit unsigned integer value
  static void bind_single_value(sqlite3_stmt* stmt, int index, std::uint32_t value) {
    sqlite3_bind_int64(stmt, index, static_cast<sqlite3_int64>(value));
  }

  // Binds a time_t value to the second placeholder in the prepared SQL
  // statement.
  // This is used for binding time values (typically representing dates or
  // timestamps).
  // stmt: Pointer to the prepared SQL statement.
  // value: time_t value to bind.
  static void bind_single_value(sqlite3_stmt* stmt, int index,
                                const std::time_t &value) {
    sqlite3_bind_int64(stmt, index, static_cast<sqlite3_int64>(value));
  }

  // Binds a double value to the third placeholder in the prepared SQL
  // statement.
  // This is typically used for binding floating-point numbers.
  // stmt: Pointer to the prepared SQL statement.
  // value: Double value to bind.
  static void bind_single_value(sqlite3_stmt* stmt, int index, double value) {
    sqlite3_bind_double(stmt, index, value);
  }

  // Binds a std::string value to the specified placeholder in the prepared SQL
  // statement.
  static void bind_single_value(sqlite3_stmt* stmt, int index,
                                const std::string &value) {
    sqlite3_bind_text(stmt, index, value.c_str(), -1, SQLITE_TRANSIENT);
  }

  // Binds a C-style string (const char*) to the specified placeholder in the
  // prepared SQL statement.
  static void bind_single_value(sqlite3_stmt* stmt, int index,
                                const char* value) {
    sqlite3_bind_text(stmt, index, value, -1, SQLITE_TRANSIENT);
  }

  template <typename T>
  static void bind_single_value(sqlite3_stmt* /*stmt*/, int /*index*/,
                                const T & /*value*/) {
    spdlog::error("Unsupported type for binding: {}", typeid(T).name());
    // exit program
    throw std::runtime_error("Unsupported type for binding");
  }

  // Template function for binding multiple values to a prepared SQL statement.
  // This function is variadic, allowing it to take multiple arguments of
  // different types. It recursively calls itself to bind each argument in the
  // parameter pack to the statement.
  // stmt: Pointer to the prepared SQL statement.
  // index: The index of the first placeholder in the statement to bind.
  // first: The first value in the parameter pack to bind.
  // rest...: The remaining values in the parameter pack.
  template <typename First, typename... Rest>
  void bind_values(sqlite3_stmt* stmt, int index, First first, Rest... rest) {
    SQLiteDatabase::bind_single_value(stmt, index, first);
    if constexpr (sizeof...(rest) > 0) {
      bind_values(stmt, index + 1, rest...);
    }
  }

public:
  // Constructor: Opens a connection to the SQLite database at the specified
  // path.
  // Throws a runtime_error if the database cannot be opened.
  explicit SQLiteDatabase(const std::string &path) {
    if (sqlite3_open_v2(path.c_str(), &db_,
                        SQLITE_OPEN_CREATE | SQLITE_OPEN_READWRITE, nullptr)
        != SQLITE_OK) {
      spdlog::error("Error opening SQLite database: {}", sqlite3_errmsg(db_));
    }
  }

  // Destructor: Closes the database connection.
  ~SQLiteDatabase() {
    if (db_ != nullptr) {
      sqlite3_close(db_);
      db_ = nullptr;
    }
  }

  // Executes a given SQL statement without expecting a return value.
  // Throws a runtime_error if the execution fails.
  void execute(const std::string &sql) {
    char* err_msg = nullptr;
    if (sqlite3_exec(db_, sql.c_str(), nullptr, nullptr, &err_msg)
        != SQLITE_OK) {
      std::string error = "SQL error msg: " + std::string(err_msg);
      sqlite3_free(err_msg);
      /* throw std::runtime_error(error); */
      spdlog::error("SQL error: {} {}", sql, error);
    }
  }

  // Prepares an SQL statement for execution and returns a pointer to the
  // prepared statement. Throws a runtime_error if the statement preparation
  // fails.
  sqlite3_stmt* prepare(const std::string &sql) {
    sqlite3_stmt* stmt = nullptr;
    if (sqlite3_prepare_v2(db_, sql.c_str(), -1, &stmt, nullptr) != SQLITE_OK) {
      /* throw std::runtime_error("Error preparing statement: " + */
      /*                          std::string(sqlite3_errmsg(db))); */
      spdlog::error("Error preparing statement: {}", sqlite3_errmsg(db_));
    }
    return stmt;
  }

  // Inserts data into the database using a prepared statement with variable
  // arguments.
  // NOTE: Curruently, insert data only support interger, time_t, and double
  //
  // Returns the ID of the inserted row
  // Throws a runtime_error on execution failure.
  template <typename... Args>
  int insert_data(const std::string &query, Args... args) {
    sqlite3_stmt* stmt = prepare(query);
    bind_values(stmt, 1, args...);

    if (sqlite3_step(stmt) != SQLITE_ROW) {
      std::string error = "Error executing insert statement: "
                          + std::string(sqlite3_errmsg(db_));
      sqlite3_finalize(stmt);
      spdlog::error(error);
      throw std::runtime_error(error);
    }

    int returned_id = sqlite3_column_int(stmt, 0);

    // Expect the next step to be SQLITE_DONE
    if (sqlite3_step(stmt) != SQLITE_DONE) {
      std::string error =
          "Execution didn't finish: " + std::string(sqlite3_errmsg(db_));
      sqlite3_finalize(stmt);
      spdlog::error(error);
      throw std::runtime_error(error);
    }

    sqlite3_finalize(stmt);
    return returned_id;
  }

  // Starts a database transaction
  void begin_transaction() { execute("BEGIN TRANSACTION;"); }

  // Commits the current transaction
  void commit_transaction() { execute("COMMIT TRANSACTION;"); }
};

class TransactionGuard {
private:
  SQLiteDatabase* db_;
  bool committed_{false};

public:
  TransactionGuard(TransactionGuard &&) = delete;
  TransactionGuard &operator=(TransactionGuard &&) = delete;
  explicit TransactionGuard(SQLiteDatabase* database) : db_(database) {
    db_->begin_transaction();
  }
  // Prevent copy and assignment
  TransactionGuard(const TransactionGuard &) = delete;
  void operator=(const TransactionGuard &) = delete;

  ~TransactionGuard() {
    if (!committed_) { db_->commit_transaction(); }
  }

  // Method to manually commit the transaction and prevent automatic commit
  void commit() {
    if (!committed_) {
      db_->commit_transaction();
      committed_ = true;
    }
  }
};

#endif
