// MMSP.data.hpp
// Class definitions for MMSP data structures
// Questions/comments to jgruber@andrew.cmu.edu (Jason Gruber)

#ifndef MMSP_DATA
#define MMSP_DATA
#include<map>
#include<cmath>
#include<vector>
#include<complex>
#include<algorithm>
#include <cstdlib>  // added ADR  9 i 09
#include <stdio.h>  // added ADR  9 i 09
#include <string.h>  // added ADR  9 i 09

namespace MMSP{

  template<typename value_type> class unused{
  public:
    // constructors
    explicit unused(...) {}

    // buffer I/O
    int buffer_size() const {return 0;}
    void to_buffer(char* buffer) const {}
    void from_buffer(const char* buffer) {}

    // assignment operators
    unused& operator=(const value_type& temp) {return *this;}
    unused& operator=(const unused& temp) {return *this;}

    // data access operators
    operator value_type&() {return static_cast<value_type>(0);}
    operator const value_type&() const {return static_cast<value_type>(0);}

    // data structure sizes
    void resize(int n) const {}
    int fields() const {return 0;}
    int nonzero() const {return 0;}
  };

  template<typename value_type> class scalar{
  public:
    // constructors
    explicit scalar(...) {}

    // buffer I/O
    int buffer_size() const {return sizeof(data);}
    void to_buffer(char* buffer) const {memcpy(buffer,&data,sizeof(data));}
    void from_buffer(const char* buffer) {memcpy(&data,buffer,sizeof(data));}

    // assignment operators
    scalar& operator=(const value_type& temp) {data=temp; return *this;}
    scalar& operator=(const scalar& temp) {data=temp.data; return *this;}

    // data access operators
    operator value_type&() {return data;}
    operator const value_type&() const {return data;}

    // data structure sizes
    void resize(int n) const {}
    int fields() const {return 1;}
    int nonzero() const {return (data!=static_cast<value_type>(0));}

  private:
    value_type data;
  };

  template<typename value_type> class vector{
  public:
    // constructors
    explicit vector(int n = 1, ...) {data.resize(n);}

    // buffer I/O
    int buffer_size() const {return data.size()*sizeof(value_type);}
    void to_buffer(char* buffer) const {memcpy(buffer,&data[0],data.size()*sizeof(value_type));}
    void from_buffer(const char* buffer) {memcpy(&data[0],buffer,data.size()*sizeof(value_type));}

    // assignment operators
    vector& operator=(const vector& temp) {data=temp.data; return *this;}

    // data access operators
    value_type& operator[](int i) {return data[i];}
    const value_type& operator[](int i) const {return data[i];}

    // data structure sizes
    void resize(int n) {return data.resize(n);}
    int fields() const {return data.size();}
    int nonzero() const {return count_if(data.begin(),data.end(),
					 bind2nd(std::not_equal_to<value_type>(),static_cast<value_type>(0)));}

  private:
    std::vector<value_type> data;
  };


  template<typename value_type> class sparse{
  public:
    // constructors
    explicit sparse(...) {}

    // buffer I/O
    int buffer_size() const;
    void to_buffer(char* buffer) const;
    void from_buffer(const char* buffer);

    // assignment operators
    sparse& operator=(const sparse& temp) {data=temp.data; return *this;}

    // data access operators
    value_type& operator[](int i);
    const value_type operator[](int i) const;

    // data structure sizes
    void resize(int n) const {}
    int fields() const {return data.size();}
    int nonzero() const {return data.size();}

    // utility functions
    int index(int i) const {return data[i].first;}
    value_type value(int i) const {return data[i].second;}

  private:
    std::vector<std::pair<int,value_type> > data;
  };

  template <typename value_type>
  value_type& sparse<value_type>::operator[](int index)
  {
    int n = data.size();
    for (int i=0; i<n; i++)
      if (data[i].first==index) return data[i].second;
    data.push_back(std::make_pair(index,static_cast<value_type>(0)));
    return data.back().second;
  }

  template <typename value_type>
  const value_type sparse<value_type>::operator[](int index) const
  {
    int n = data.size();
    for (int i=0; i<n; i++)
      if (data[i].first==index) return data[i].second;
    return static_cast<value_type>(0);
  }

  template <typename value_type>
  int sparse<value_type>::buffer_size() const
  {
    int n = data.size();
    return sizeof(int)+sizeof(std::pair<int,value_type>);
  }

  template <typename value_type>
  void sparse<value_type>::to_buffer(char* buffer) const
  {
    char* p = buffer;
    int n = data.size();
    memcpy(p,&n,sizeof(int));
    p += sizeof(int);
    memcpy(p,&data[0],n*sizeof(std::pair<int,value_type>));
  }

  template <typename value_type>
  void sparse<value_type>::from_buffer(const char* buffer) 
  {
    const char* p = buffer;
    int n;
    memcpy(&n,p,sizeof(int));
    p += sizeof(int);
    data.resize(n);
    memcpy(&data[0],p,n*sizeof(std::pair<int,value_type>));
  }

} // namespace MMSP

#endif
