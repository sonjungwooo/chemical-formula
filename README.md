# chemical_formula 
Made for sci-competition

## 탐구 주제
프로그래밍을 통하여 우리가 평소에 외우기 힘들었던 화학식을 자동으로 완성해주어 더 쉽고 빠르게 암기를 할 수 있도록 하는 앱을 만드는 것

##탐구하게 된 동기
저희 또래의 아이들은 중학교 2학년에 분자와 원자를 배우고, 3학년 때 본격적으로 화학반응식을 외우기 시작합니다.
하지만 이런 화학반응식을 외우는데 어려움을 느끼는 친구를 자주 보았고 또한 저희도 처음에 외울 때는 헷갈렸습니다.
하지만 이 프로그램을 만듦으로써 화학반응식을 좀더 쉽고 편리하게 외울 수 있게 하려고 합니다.

##선행 연구 및 고찰 
화학 반응식은 화학 반응이 일어나면 반응하는 물질인 반응물과 반응물 사이에서 원자의 재배열이 일어나 새로운 결합을 형성하게 된 생성물을 화학식, 기호를 이용하
여 나타낸 식입니다.
이 결합과정 중 하나가 이온결합 인데, 이것이 우리가 학교에서 배우는 것입니다. 하지만 이 화학반응식을 쓸 때 각 원소의  전하량이 같지 않으므로 원자의 개수를 
다르게 해주어야 하고 , 그 과정에서 각 분자들의 계수를 맞춰줘야 합니다. 
중학교 교과서에는 이온결합 밖에 실려 있지 않지만, 이 이외에도 금속결합, 공유결합이 있습니다.

##탐구를 실행한 절차 
어떤 프로그램을 사용해 만들 것인지를 생각합니다.
어떤 프로그래밍 언어를 사용해 만들 것인지를 생각합니다,
먼저 프로그램을 짜기 전에 어떻게 프로그램을 짤지 구상합니다:
사용자가 입력한 화학반응식을 받는 프로그램을 짜고, 그 받은 화학식을 어떻게 나눌지도 구상합니다,
그리고 그 나눈 화학식의 원자들을 각각 분석해 그에 맞는 분자의 계수비를 만들어 출력하도록 구상합니다.
프로그램 구상이 끝나면 프로그램을 만듭니다.
만드는 도중에 오류나 에러가 뜨면 다시 시행착오를 겪으면서 오류를 고쳐 나갑니다,
오류 없이 전부 작동되는 것을 확인합니다,
프로그램을 여러 번 실행해보면서 오류를 다시 찾습니다.

##탐구 방법
프로그래밍 언어는 C++이라는 언어를 사용하였고 비주얼 스튜디오라는 앱을 통해서 프로그래밍을 했습니다.
비주얼 스튜디오에서 프로그래밍을 하고 수시로 실행을 하여 오류가 있는지 확인하고 화학반응식에서의 계수비의를 이용해 연립방적식을 만들어 프로그램을
만듭니다.
이 만든 프로그램을 다시 실행해보는 단계에서 많은 경우의 수가 나올 수 있는 화학 반응식을 입력함으로써 오류가 있는지 없는지 찾아낼 수 있고, 많이 시도하여
오류를 최대한 줄여가며 프로그래밍 했습니다.
