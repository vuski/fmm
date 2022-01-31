# Fast Map Matching

- 이 repository는 Can Yang의 https://github.com/cyang-kth/fmm 프로젝트에서 일부 내용을 수정한 것입니다.
- GDAL이나 Boost등의 라이브러리를 덜어내고 header-only 로 정리했습니다.
- 물론 그 과정에서 shape 파일을 읽는 부분 등은 삭제되었습니다. geojson 을 읽습니다.
- 실행에 필요한 third-party 라이브러리들도 모두 include 폴더에 있습니다. 
- 윈도우의 visual studio 에서만 테스트 되었습니다. sln 파일을 열면 로딩되는 프로젝트는 release - x64 모드에서 실행됩니다.
- 기타 자세한 내용은 https://www.vw-lab.com/104   에 설명되어 있습니다.
- 원 저작자의 리포지토리와 마찬가지로 apache license 2.0 을 따릅니다.