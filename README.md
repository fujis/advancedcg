# アドバンストCG / Advanced CG
筑波大学情報学群情報メディア創成学類・情報科学類の授業「アドバンストCG」の第5回～第8回レポートのためのベースとなるコード公開用リポジトリです．
第1回～第4回については金森先生の[リポジトリ](https://github.com/yshhrknmr/AdvancedCG)を参照してください．
This is a git repository for the lecture "Advanced Computer Graphics" at the University of Tsukuba, School of Informatics, College of Media Arts, Science and Technology. It includes the base code written in C++ for the reports.

## フォルダ構成など
- Advanced05～Advanced08: 第5回～第8回の授業課題のベースプログラムが置いてあります．中にはソースコードとVisual Studio用のプロジェクトファイルなどが入っています．Visual Studioから開く場合はこのフォルダ内のプロジェクトファイルではなく，リポジトリ直下に置いてある Advanced05Assignment_vs2026.sln (Visual Studio 2026用) または Advanced05Assignment_vs2022.sln (Visual Studio 2022用) を開いてください．
- bin: 実行時フォルダで，実行時に参照するモデルファイルや画像ファイルなどもここに置いてあります．
- include: 第5回～第8回の授業課題で必要となるヘッダファイルが置いてあります．Visual Studio 以外でコンパイルする際はこのディレクトリ (include) をインクルードフォルダとして追加してください．
- vc-x64-lib-common: Visual Studio 共通のライブラリファイル
- vc-x64-libs: Visual Studio のバージョン別のライブラリファイル
- licenses: 使用している外部ライブラリのライセンス文

## 使用している外部ライブラリについて
各ライブラリのライセンス文は licenses フォルダにまとめてあります(stb_image, stb_image_write はヘッダファイルの末尾に記載)．
- [GLEW](https://github.com/nigels-com/glew) (include/GL, vc-x64-lib-common/glew32.*, bin/glew32.dll): OpenGL拡張機能の読み込み．Modified BSD License (一部 Mesa の MIT License，Khronos の MIT License)．→ licenses/GLEW_LICENSE.txt
- [GLFW](https://www.glfw.org/) (include/GLFW, vc-x64-libs, bin/glfw3.dll): ウィンドウ作成と入力処理．zlib/libpng License．→ licenses/GLFW_LICENSE.txt
- [GLM](https://github.com/g-truc/glm) 0.9.4 (include/glm): ベクトル・行列計算．MIT License．→ licenses/GLM_LICENSE.txt
- [Dear ImGui](https://github.com/ocornut/imgui) 1.83 (Advanced05～08/imgui): GUI．MIT License．→ licenses/IMGUI_LICENSE.txt
- [stb_image / stb_image_write](https://github.com/nothings/stb) (include/stb_image.h, include/stb_image_write.h): 画像ファイルの読み書きに使用(第5回)．public domain (Unlicense) と MIT License のどちらかを選べる二重ライセンスで，ライセンス文は各ヘッダファイルの末尾に記載されています．

## トラブルシューティング
 - Visual Studioはgithub連携機能がありますが，このリポジトリには複数のVisual Studioソリューションとソースコード以外のファイルが多く含まれているのでうまくいかない可能性があります．別途リポジトリをクローンするかDownload ZIPで全体をダウンロード・解凍後，リポジトリ直下の*.slnファイルをVisual Studioで開くようにしてください．

 - Visual Studioで「error MSB8036: Windows SDK バージョン 10.0.xxxxxx.0 が見つかりませんでした。」というエラーが出てビルドできない．  
　⇒ 「プロジェクト」メニューから「プロジェクトの再ターゲット」で「Windows SDK バージョン:」 のところにその環境で対応するバージョンが出るので，問題なければそのまま「OK」をクリック
 - Visual Studio 2026で_vs2022.slnを開くとプラットフォームツールセットの再ターゲットを求められることがあります．その場合は_vs2026.slnの方を開いてください．
 - Visual Studioで実行しても特定のプロジェクトしか実行されない場合は，「プロジェクト」メニューから「スタートアップ プロジェクトの構成」で「現在の選択」にチェックしてください．

 - 自身のPC環境(Windows)で動かしたい人向けに講義資料フォルダ(SharePointの方)に「Visual Studio Community インストール方法」というPDFを用意してあるので参考にしてください．
