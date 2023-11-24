# advancedcg_report
筑波大学情報学群情報メディア創成学類の授業「アドバンストCG」の第5回～第8回のレポートのためのベースコード公開用リポジトリです．
This is a repository for the lecture course "Advanced Computer Graphics" in University of Tsukuba. 

## 各フォルダについて
* Advanced05～Advanced08 : 第5回～第8回の授業課題のベースプログラムが置いてあります．
* bin : 実行用フォルダ．DLLや実行時に参照する画像やモデルファイル，シェーダーなどが置いてあります．
* include : 第5回～第8回の授業課題で必要となるヘッダファイルが置いてあります．VisualStudio以外でコンパイルする際はincludeをインクルードフォルダとして追加するか，include内の*.hファイルをフォルダ構造そのままにcppファイルのあるところにコピーするようにしてください．
* vc-x64-lib-common : Visual Studio 共通のライブラリファイル
* vc-x64-libs : バージョン別のライブラリファイル

## 注意事項
### Windows
* もし自分の PC に Visual Studio (2019 または 2022) がインストールされていない場合は，講義資料フォルダ(SharePoint)の「Visual Studio Community インストール方法.pdf」 を参考にしてインストールしてください．
* 上記ソースコード一式に含まれるVisualStudioソリューション/プロジェクトファイルを用いた場合で，「error MSB8036: Windows SDK バージョン 10.0.xxxxxx.0 が見つかりませんでした。」というエラーが出てビルドできない場合は，「プロジェクト」メニューから「プロジェクトの再ターゲット」を実行してください．「Windows SDK バージョン:」 のところにその環境で対応するバージョンが出るので，問題なければそのまま「OK」をクリックすればプロジェクトがその環境に合わせて更新されます．
* 全学計算機システムのリモートデスクトップでプログラムを実行した場合，第5回の課題用プログラムでメッシュおよびテクスチャがうまく表示されないようです(実習室のPCでは問題なし)．これは環境が対応しているGLSLのバージョン違いが原因です．リモートデスクトップで第5回の課題を行う場合は， bin/shaders/shading.fp と bin/shaders/shading.vp の12or13行目の
`#version 410 core`
というところを
`#version 330`
と書き換えてください．

### Mac, Linux (Ubuntu)
* 上記ソースコード一式はMac用MakefileとUbuntu用Makefileは含んでいるのでそれを使ってください(Ubuntu用はMakefile.ubuntuとしてあるので適宜ファイル名変更してください)．ただし，画像などのファイルはincludeフォルダなどと同じレベルのbinフォルダを作成して，その中に入れるようにしてください．
* 第5-8回のMakefileは全て共通です．実行ファイル名がすべて"acg"と設定されているので注意してください．課題毎に変えたい場合はMakefileの13行目:
TARGETS = acg
のacgのところを随時変更してください．また，自身のMacを使って課題を行う場合で，M1チップ搭載Macの場合でHomebrewを使っている場合は，6行目(LDFLAGS)で -L/opt/local/lib としているところを -L/opt/homebrew/lib に変更してみてください．
* Mac環境でマウスによる制御点選択の位置がずれる可能性があります．その場合，アプリケーション起動後にウィンドウサイズを変えると直るようです．
